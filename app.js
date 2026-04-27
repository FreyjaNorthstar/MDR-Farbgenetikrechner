const LOCI = [
  // Reihenfolge passend zur UI:
  // Extension, Agouti, Dun, (Cream/Pearl), Champagne, Grey, Silver
  {
    key: "E",
    label: "Extension (E/e)",
    alleleOrder: ["E", "e"],
    genotypes: ["EE", "Ee", "ee"],
  },
  {
    key: "A",
    label: "Agouti (A+, A, At, a)",
    alleleOrder: ["A+", "A", "At", "a"], // Dominanz: A+ > A > At > a
    genotypes: ["A+A+", "A+A", "A+At", "A+a", "AA", "AAt", "Aa", "AtAt", "Ata", "aa"],
  },
  {
    key: "D",
    label: "Dun (D/d)",
    alleleOrder: ["D", "d"],
    genotypes: ["DD", "Dd", "dd"],
  },
  {
    key: "CrPrl",
    label: "Cream / Pearl (Cr/cr/prl)",
    alleleOrder: ["Cr", "cr", "prl"],
    // mögliche Kombis: crcr, crprl, Crcr, Crprl, CrCr, prlprl
    // (Reihenfolge wird intern via alleleOrder normalisiert)
    genotypes: ["CrCr", "Crcr", "Crprl", "crcr", "crprl", "prlprl"],
  },
  {
    key: "Ch",
    label: "Champagne (Ch/ch)",
    alleleOrder: ["Ch", "ch"],
    genotypes: ["ChCh", "Chch", "chch"],
  },
  {
    key: "G",
    label: "Grey (G/g)",
    alleleOrder: ["G", "g"],
    genotypes: ["GG", "Gg", "gg"],
  },
  {
    key: "Z",
    label: "Silver (Z/z)",
    alleleOrder: ["Z", "z"],
    genotypes: ["ZZ", "Zz", "zz"],
  },
];

function $(id) {
  const el = document.getElementById(id);
  if (!el) throw new Error(`Element nicht gefunden: ${id}`);
  return el;
}

function normalizeGenotype(genotype, alleleOrder) {
  // genotype: z.B. "eE" -> "Ee", "crCr" -> "Crcr"
  // alleleOrder: ["E","e"] oder ["Cr","cr"]
  const alleles = splitAlleles(genotype, alleleOrder);
  const sorted = alleles.sort((a, b) => alleleOrder.indexOf(a) - alleleOrder.indexOf(b));
  return sorted.join("");
}

function splitAlleles(genotype, alleleOrder) {
  // Unterstützt 1- oder 2-Zeichen-Allele (z.B. "Cr" vs "cr").
  // Wir parsen, indem wir mit den bekannten Allelen matchen.
  const candidates = [...alleleOrder].sort((a, b) => b.length - a.length);
  const alleles = [];
  let rest = genotype;
  while (rest.length > 0) {
    const match = candidates.find((a) => rest.startsWith(a));
    if (!match) {
      throw new Error(`Genotyp konnte nicht geparst werden: "${genotype}" (Rest: "${rest}")`);
    }
    alleles.push(match);
    rest = rest.slice(match.length);
  }
  if (alleles.length !== 2) {
    throw new Error(`Genotyp muss aus 2 Allelen bestehen: "${genotype}"`);
  }
  return alleles;
}

function gametesFromGenotype(genotype, alleleOrder) {
  const [a1, a2] = splitAlleles(genotype, alleleOrder);
  if (a1 === a2) return [{ allele: a1, p: 1 }];
  return [
    { allele: a1, p: 0.5 },
    { allele: a2, p: 0.5 },
  ];
}

function punnett(parent1Genotype, parent2Genotype, alleleOrder) {
  const g1 = normalizeGenotype(parent1Genotype, alleleOrder);
  const g2 = normalizeGenotype(parent2Genotype, alleleOrder);
  const gam1 = gametesFromGenotype(g1, alleleOrder);
  const gam2 = gametesFromGenotype(g2, alleleOrder);

  const dist = new Map(); // genotype -> probability
  for (const a of gam1) {
    for (const b of gam2) {
      const child = normalizeGenotype(`${a.allele}${b.allele}`, alleleOrder);
      const p = a.p * b.p;
      dist.set(child, (dist.get(child) ?? 0) + p);
    }
  }
  return dist;
}

function cartesianCombineDistributions(distsByLocus) {
  // Input: [{locusKey, dist: Map(genotype->p)}]
  // Output: Map(signature->p) where signature like "EE|Aa|gg|Crcr|Dd" in locus order
  let acc = new Map([["", 1]]);

  for (const { locusKey, dist } of distsByLocus) {
    const next = new Map();
    for (const [sig, pSig] of acc.entries()) {
      for (const [gt, pGt] of dist.entries()) {
        const newSig = sig ? `${sig}|${gt}` : gt;
        next.set(newSig, (next.get(newSig) ?? 0) + pSig * pGt);
      }
    }
    acc = next;
  }
  return acc;
}

function roundPct(x) {
  return Math.round(x * 1000) / 10; // 0.1%
}

function agoutiCategory(agoutiGenotype) {
  // Rückgabe: "wildbay" | "bay" | "sealbrown" | "black"
  // Dominanz: A+ > A > At > a
  if (agoutiGenotype === "aa") return "black";
  if (agoutiGenotype.includes("A+")) return "wildbay";
  // Wichtig: "A+" enthält auch "A" als Teilstring, daher A+ zuerst prüfen.
  if (agoutiGenotype.includes("A")) return "bay";
  if (agoutiGenotype.includes("At")) return "sealbrown";
  return "black";
}

function derivePhenotype({ E, A, G, CrPrl, Ch, Z, D }) {
  // Sehr vereinfachte Ableitung, aber konsistent und erweiterbar.
  // 1) Basisfarbe via E/A
  const hasBlackPigment = E !== "ee";
  let base;
  if (!hasBlackPigment) {
    base = "Chestnut";
  } else {
    const aCat = agoutiCategory(A);
    if (aCat === "black") base = "Black";
    else if (aCat === "wildbay") base = "Wildbay";
    else if (aCat === "bay") base = "Bay";
    else base = "Sealbrown";
  }

  // 2) Cream / Pearl (gleicher Locus)
  const crCount = CrPrl === "CrCr" ? 2 : CrPrl === "Crcr" || CrPrl === "Crprl" ? 1 : 0;
  const prlCount = CrPrl === "prlprl" ? 2 : CrPrl === "crprl" || CrPrl === "Crprl" ? 1 : 0;

  const isPearl = prlCount === 2;
  const isPearlCarrier = prlCount === 1 && crCount === 0; // crprl
  const isCreamPearl = prlCount === 1 && crCount === 1; // Crprl

  let creamNote = null;
  let pearlNote = null;
  let baseLight0 = base;

  // Cream-Benennung (nur wenn nicht prlprl und nicht Crprl)
  if (!isPearl && !isCreamPearl) {
    if (crCount === 1) {
      if (base === "Chestnut") baseLight0 = "Palomino";
      else if (base === "Bay" || base === "Wildbay") baseLight0 = "Buckskin";
      else if (base === "Black") baseLight0 = "Smoky Black";
      else if (base === "Sealbrown") baseLight0 = "Smoky Brown";
      creamNote = "Cream (Crcr)";
    } else if (crCount === 2) {
      if (base === "Chestnut") baseLight0 = "Cremello";
      else if (base === "Bay" || base === "Wildbay") baseLight0 = "Perlino";
      else if (base === "Black") baseLight0 = "Smoky Cream";
      else if (base === "Sealbrown") baseLight0 = "Sealbrown Cream";
      creamNote = "Cream (CrCr)";
    } else {
      baseLight0 = base;
    }
  }

  // Pearl (prlprl)
  if (isPearl) {
    if (base === "Chestnut") baseLight0 = "Apricot";
    else if (base === "Bay" || base === "Wildbay") baseLight0 = "Pearl Bay";
    else if (base === "Black") baseLight0 = "Pearl Black";
    else if (base === "Sealbrown") baseLight0 = "Pearl Brown";
    else baseLight0 = `${base} (Pearl)`;
    pearlNote = "Pearl (prlprl)";
  } else if (isCreamPearl) {
    // Crprl wird wie CrCr benannt
    if (base === "Chestnut") baseLight0 = "Cremello";
    else if (base === "Bay" || base === "Wildbay") baseLight0 = "Perlino";
    else if (base === "Black") baseLight0 = "Smoky Cream";
    else if (base === "Sealbrown") baseLight0 = "Sealbrown Cream";
    else baseLight0 = `${base} (Crprl)`;
    // im Namen nicht extra berücksichtigen, aber als Tag ist es manchmal hilfreich
    pearlNote = "Cream/Pearl (Crprl)";
  }

  // 3) Champagne (vereinfacht)
  const chCount = Ch === "ChCh" ? 2 : Ch === "Chch" ? 1 : 0;
  let champagneNote = null;
  let baseLight = baseLight0;
  if (chCount > 0) {
    champagneNote = chCount === 2 ? "Champagne (homozygot)" : "Champagne";
    // Champagne-Benennung basiert auf der Grundfarbe (E/A)
    if (base === "Chestnut") baseLight = "Gold Champagne";
    else if (base === "Bay" || base === "Wildbay") baseLight = "Amber Champagne";
    else if (base === "Black") baseLight = "Classic Champagne";
    else if (base === "Sealbrown") baseLight = "Sable Champagne";
    else baseLight = `${baseLight0} (Champagne)`;
  }

  // 4) Silver (vereinfacht)
  const zCount = Z === "ZZ" ? 2 : Z === "Zz" ? 1 : 0;
  const isSilver = zCount > 0 && hasBlackPigment;
  let silverNote = null;
  let baseSilver = baseLight;
  if (zCount > 0) {
    if (hasBlackPigment) {
      silverNote = zCount === 2 ? "Silver (homozygot)" : "Silver";
      baseSilver = `${baseLight} (Silver)`;
    } else {
      silverNote = "Silver (ohne Wirkung bei Fuchs)";
    }
  }

  // 5) Dun
  const isDun = D !== "dd";
  let dunNote = null;
  let withDun = baseSilver;
  if (isDun) {
    dunNote = "Dun";
    // gewünschte Benennung
    if (base === "Chestnut") withDun = "Red Dun";
    else if (base === "Bay" || base === "Wildbay") withDun = "Classic Dun";
    else if (base === "Black") withDun = "Grulla";
    else if (base === "Sealbrown") withDun = "Brown Dun";
    else withDun = `${baseSilver} Dun`;
  }

  // 6) Grey überschreibt (langfristig)
  const isGrey = G !== "gg";
  const tags = [];
  if (creamNote) tags.push(creamNote);
  if (pearlNote) tags.push(pearlNote);
  else if (isPearlCarrier) tags.push("Pearl Träger (crprl)");
  if (champagneNote) tags.push(champagneNote);
  if (silverNote) tags.push(silverNote);
  if (dunNote) tags.push(dunNote);
  if (isGrey) tags.push("Grey");

  const shown = isGrey ? "Grey" : withDun;
  const detail = isGrey ? `Grundfarbe darunter: ${withDun}` : null;
  return { shown, detail, tags };
}

function formatDistributionTable(entries, caption) {
  const rows = entries
    .sort((a, b) => b.p - a.p)
    .map(
      (e) =>
        `<tr><td>${escapeHtml(e.label)}</td><td>${roundPct(e.p)}%</td></tr>`
    )
    .join("");
  return `
    <div class="pill">${escapeHtml(caption)}</div>
    <table>
      <thead><tr><th>Ausprägung</th><th>Wahrscheinlichkeit</th></tr></thead>
      <tbody>${rows}</tbody>
    </table>
  `;
}

function escapeHtml(s) {
  return String(s)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function readParents() {
  const p1 = {};
  const p2 = {};

  for (const locus of LOCI) {
    p1[locus.key] = $(makeSelectId(1, locus.key)).value;
    p2[locus.key] = $(makeSelectId(2, locus.key)).value;
  }
  return { p1, p2 };
}

function makeSelectId(parentIdx, locusKey) {
  return `p${parentIdx}_${locusKey}`;
}

function render() {
  const { p1, p2 } = readParents();
  const distsByLocus = LOCI.map((locus) => ({
    locusKey: locus.key,
    dist: punnett(p1[locus.key], p2[locus.key], locus.alleleOrder),
  }));

  const combined = cartesianCombineDistributions(distsByLocus);
  const locusKeys = LOCI.map((l) => l.key);

  // Phänotypen zusammenfassen
  const phenoDist = new Map(); // shown -> p
  const detailDist = new Map(); // detail string -> p (optional)

  for (const [sig, p] of combined.entries()) {
    const parts = sig.split("|");
    const g = {};
    for (let i = 0; i < locusKeys.length; i++) g[locusKeys[i]] = parts[i];

    const ph = derivePhenotype(g);
    phenoDist.set(ph.shown, (phenoDist.get(ph.shown) ?? 0) + p);
    if (ph.detail) detailDist.set(ph.detail, (detailDist.get(ph.detail) ?? 0) + p);
  }

  // Einzel-Locus Tabellen
  const locusTables = distsByLocus
    .map(({ locusKey, dist }) => {
      const locus = LOCI.find((l) => l.key === locusKey);
      const entries = [...dist.entries()].map(([gt, p]) => ({ label: gt, p }));
      return formatDistributionTable(entries, `${locusKey}: ${locus?.label ?? locusKey}`);
    })
    .join("");

  const phenoEntries = [...phenoDist.entries()].map(([label, p]) => ({ label, p }));
  const phenoTable = formatDistributionTable(phenoEntries, "Phänotyp (vereinfacht)");

  const details = [...detailDist.entries()]
    .sort((a, b) => b[1] - a[1])
    .slice(0, 6)
    .map(([d, p]) => `<li>${escapeHtml(d)}: ${roundPct(p)}%</li>`)
    .join("");

  const detailsBlock =
    details.length > 0
      ? `<div class="muted" style="margin-top:10px">Hinweis (nur bei Grey):</div><ul>${details}</ul>`
      : "";

  $("result").classList.remove("muted");
  $("result").innerHTML = `
    ${phenoTable}
    ${detailsBlock}
    <div style="height:12px"></div>
    <div class="pill">Genotyp-Verteilungen je Locus</div>
    ${locusTables}
  `;
}

function resetToDefaults() {
  // sinnvolle Defaultwerte: heterozygot bei Dominanzloci, damit man direkt Variation sieht
  const defaults = {
    E: "Ee",
    A: "Aa",
    G: "gg",
    CrPrl: "crcr",
    Ch: "chch",
    Z: "zz",
    D: "dd",
  };

  for (const locus of LOCI) {
    $(makeSelectId(1, locus.key)).value = defaults[locus.key] ?? locus.genotypes[0];
    $(makeSelectId(2, locus.key)).value = defaults[locus.key] ?? locus.genotypes[0];
  }
  $("result").classList.add("muted");
  $("result").textContent = "Wähle die Genotypen und klicke auf „Berechnen“.";
}

function init() {
  for (const locus of LOCI) {
    for (const parentIdx of [1, 2]) {
      const sel = $(makeSelectId(parentIdx, locus.key));
      sel.innerHTML = locus.genotypes
        .map((gt) => `<option value="${gt}">${gt}</option>`)
        .join("");
    }
  }

  $("calcBtn").addEventListener("click", render);
  $("resetBtn").addEventListener("click", resetToDefaults);

  resetToDefaults();
}

document.addEventListener("DOMContentLoaded", init);

