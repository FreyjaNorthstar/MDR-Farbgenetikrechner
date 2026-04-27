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
    key: "Cr",
    label: "Cream (Cr/cr)",
    alleleOrder: ["Cr", "cr"],
    genotypes: ["CrCr", "Crcr", "crcr"],
  },
  {
    key: "Prl",
    label: "Pearl (Prl/prl)",
    alleleOrder: ["Prl", "prl"],
    genotypes: ["PrlPrl", "Prlprl", "prlprl"],
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

function derivePhenotype({ E, A, G, Cr, Prl, Ch, Z, D }) {
  // Sehr vereinfachte Ableitung, aber konsistent und erweiterbar.
  // 1) Basisfarbe via E/A
  const hasBlackPigment = E !== "ee";
  let base;
  if (!hasBlackPigment) {
    base = "Fuchs";
  } else {
    const aCat = agoutiCategory(A);
    if (aCat === "black") base = "Rappe";
    else if (aCat === "wildbay") base = "Wildbay";
    else if (aCat === "bay") base = "Brauner";
    else base = "Schwarzbraun (Seal Brown)";
  }

  // 2) Cream
  const crCount = Cr === "CrCr" ? 2 : Cr === "Crcr" ? 1 : 0;
  let creamNote = null;
  let baseCream = base;
  if (crCount === 1) {
    if (base === "Fuchs") baseCream = "Palomino";
    if (base === "Brauner") baseCream = "Buckskin";
    if (base === "Rappe") baseCream = "Smoky Black";
    creamNote = "1× Cream";
  } else if (crCount === 2) {
    if (base === "Fuchs") baseCream = "Cremello";
    if (base === "Brauner") baseCream = "Perlino";
    if (base === "Rappe") baseCream = "Smoky Cream";
    creamNote = "2× Cream";
  }

  // 3) Pearl (vereinfacht) + Interaktion mit Cream
  const prlCount = Prl === "PrlPrl" ? 2 : Prl === "Prlprl" ? 1 : 0;
  const isPearl = prlCount === 2;
  const isPearlCarrier = prlCount === 1;

  let pearlNote = null;
  let basePearl = baseCream;

  if (isPearl) {
    pearlNote = "Pearl";

    // Cream+Pearl wird oft als "pseudo double cream" beschrieben.
    // Für dieses MVP modellieren wir das so: 1×Cream + prlprl => wie 2×Cream in der Benennung.
    if (crCount === 1) {
      if (base === "Fuchs") basePearl = "Cremello (Cream+Pearl)";
      else if (base === "Brauner") basePearl = "Perlino (Cream+Pearl)";
      else if (base === "Rappe") basePearl = "Smoky Cream (Cream+Pearl)";
      else basePearl = `${baseCream} (Cream+Pearl)`;
      pearlNote = "Pearl (mit Cream)";
    } else if (crCount === 0) {
      // Ohne Cream benennen wir Pearl separat, angelehnt an gängige Begriffe.
      if (base === "Fuchs") basePearl = "Apricot Pearl";
      else if (base === "Brauner") basePearl = "Amber Pearl";
      else if (base === "Rappe") basePearl = "Sable Pearl";
      else basePearl = `${baseCream} (Pearl)`;
    } else {
      // crCount === 2: Name bleibt, Pearl nur als Tag
      basePearl = baseCream;
      pearlNote = "Pearl (mit 2×Cream)";
    }
  }

  // 4) Champagne (vereinfacht)
  const chCount = Ch === "ChCh" ? 2 : Ch === "Chch" ? 1 : 0;
  let champagneNote = null;
  let baseLight = basePearl;
  if (chCount > 0) {
    champagneNote = chCount === 2 ? "Champagne (homozygot)" : "Champagne";
    // Benennung basiert hier auf der *Grundfarbe* (E/A), weil das am stabilsten ist.
    // Kombis (z.B. Cream+Champagne) werden als Tag geführt und nicht als eigene Namensmatrix.
    if (base === "Fuchs") baseLight = "Gold Champagne";
    else if (base === "Brauner") baseLight = "Amber Champagne";
    else if (base === "Rappe") baseLight = "Classic Champagne";
    else baseLight = `${basePearl} (Champagne)`;
  }

  // 5) Silver (vereinfacht)
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

  // 6) Dun
  const isDun = D !== "dd";
  let dunNote = null;
  let withDun = baseSilver;
  if (isDun) {
    dunNote = "Dun";
    // grobe deutsche Bezeichnung
    if (baseSilver === "Fuchs") withDun = "Fuchsfalbe (Red Dun)";
    else if (baseSilver === "Brauner") withDun = "Falbe (Bay Dun)";
    else if (baseSilver === "Rappe") withDun = "Mausfalbe (Grullo)";
    else withDun = `${baseSilver} (Dun)`;
  }

  // 7) Grey überschreibt (langfristig)
  const isGrey = G !== "gg";
  const tags = [];
  if (creamNote) tags.push(creamNote);
  if (pearlNote) tags.push(pearlNote);
  else if (isPearlCarrier) tags.push("Pearl (Träger)");
  if (champagneNote) tags.push(champagneNote);
  if (silverNote) tags.push(silverNote);
  if (dunNote) tags.push(dunNote);
  if (isGrey) tags.push("Schimmel (Grey)");

  const shown = isGrey ? "Schimmel (Grey)" : withDun;
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

function parseCreamPearlCombo(value) {
  // value: "Crcr|prlprl" => { Cr: "Crcr", Prl: "prlprl" }
  const [Cr, Prl] = String(value).split("|");
  if (!Cr || !Prl) throw new Error(`Ungültige Cream/Pearl-Kombination: "${value}"`);
  return { Cr, Prl };
}

function readParents() {
  const p1 = {};
  const p2 = {};

  // Cream/Pearl kommen aus einem gemeinsamen Dropdown
  const c1 = parseCreamPearlCombo($("p1_CrPrl").value);
  const c2 = parseCreamPearlCombo($("p2_CrPrl").value);
  p1.Cr = c1.Cr;
  p1.Prl = c1.Prl;
  p2.Cr = c2.Cr;
  p2.Prl = c2.Prl;

  for (const locus of LOCI) {
    // Cr/Prl sind bereits gesetzt (gemeinsames Dropdown)
    if (locus.key === "Cr" || locus.key === "Prl") continue;
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
    Cr: "crcr",
    Prl: "prlprl",
    Ch: "chch",
    Z: "zz",
    D: "dd",
  };

  // gemeinsames Cream/Pearl Default
  $("p1_CrPrl").value = `${defaults.Cr}|${defaults.Prl}`;
  $("p2_CrPrl").value = `${defaults.Cr}|${defaults.Prl}`;

  for (const locus of LOCI) {
    if (locus.key === "Cr" || locus.key === "Prl") continue;
    $(makeSelectId(1, locus.key)).value = defaults[locus.key] ?? locus.genotypes[0];
    $(makeSelectId(2, locus.key)).value = defaults[locus.key] ?? locus.genotypes[0];
  }
  $("result").classList.add("muted");
  $("result").textContent = "Wähle die Genotypen und klicke auf „Berechnen“.";
}

function init() {
  for (const locus of LOCI) {
    for (const parentIdx of [1, 2]) {
      if (locus.key === "Cr" || locus.key === "Prl") continue;
      const sel = $(makeSelectId(parentIdx, locus.key));
      sel.innerHTML = locus.genotypes
        .map((gt) => `<option value="${gt}">${gt}</option>`)
        .join("");
    }
  }

  // gemeinsames Cream/Pearl Dropdown füllen (3×3 Kombis)
  const cr = LOCI.find((l) => l.key === "Cr")?.genotypes ?? ["CrCr", "Crcr", "crcr"];
  const prl = LOCI.find((l) => l.key === "Prl")?.genotypes ?? ["PrlPrl", "Prlprl", "prlprl"];
  const combos = [];
  for (const c of cr) for (const p of prl) combos.push({ value: `${c}|${p}`, label: `${c} + ${p}` });
  const comboHtml = combos.map((o) => `<option value="${o.value}">${o.label}</option>`).join("");
  $("p1_CrPrl").innerHTML = comboHtml;
  $("p2_CrPrl").innerHTML = comboHtml;

  $("calcBtn").addEventListener("click", render);
  $("resetBtn").addEventListener("click", resetToDefaults);

  resetToDefaults();
}

document.addEventListener("DOMContentLoaded", init);

