const LOCI = [
  // Reihenfolge passend zur UI:
  // Extension, Agouti, Dun, (Cream/Pearl), Champagne, Grey, Silver, Overo, Splashed White, Flaxen, Sooty, Rabicano
  {
    key: "E",
    label: "Extension (E/e)",
    alleleOrder: ["E", "e"],
    genotypes: ["EE", "Ee", "ee"],
  },
  {
    key: "A",
    label: "Agouti (Ap, A1, At, a0)",
    alleleOrder: ["Ap", "A1", "At", "a0"], // Dominanz: Ap > A1 > At > a0
    genotypes: [
      "ApAp",
      "ApA1",
      "ApAt",
      "Apa0",
      "A1A1",
      "A1At",
      "A1a0",
      "AtAt",
      "Ata0",
      "a0a0",
    ],
  },
  {
    key: "D",
    label: "Dun (D/d)",
    alleleOrder: ["D", "d"],
    genotypes: ["DD", "Dd", "dd"],
  },
  {
    key: "CrPrl",
    label: "Cream / Pearl (Cr/cr/pl)",
    alleleOrder: ["Cr", "cr", "pl"],
    // mögliche Kombis: crcr, crprl, Crcr, Crprl, CrCr, prlprl
    // (Reihenfolge wird intern via alleleOrder normalisiert)
    genotypes: ["CrCr", "Crcr", "Crpl", "crcr", "crpl", "plpl"],
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
  {
    key: "O",
    label: "Overo (O/o)",
    alleleOrder: ["O", "o"],
    genotypes: ["OO", "Oo", "oo"],
  },
  {
    key: "SPL",
    label: "Splashed White (SPL/spl)",
    alleleOrder: ["SPL", "spl"],
    genotypes: ["SPLSPL", "SPLspl", "splspl"],
  },
  {
    key: "Fl",
    label: "Flaxen (Fl/fl)",
    alleleOrder: ["Fl", "fl"],
    genotypes: ["FlFl", "Flfl", "flfl"],
  },
  {
    key: "Sty",
    label: "Sooty (Sty/sty)",
    alleleOrder: ["Sty", "sty"],
    genotypes: ["StySty", "Stysty", "stysty"],
  },
  {
    key: "Ra",
    label: "Rabicano (Ra/ra)",
    alleleOrder: ["Ra", "ra"],
    genotypes: ["RaRa", "Rara", "rara"],
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
  // Dominanz: Ap > A1 > At > a0
  if (agoutiGenotype === "a0a0") return "black";
  if (agoutiGenotype.includes("Ap")) return "wildbay";
  if (agoutiGenotype.includes("A1")) return "bay";
  if (agoutiGenotype.includes("At")) return "sealbrown";
  return "black";
}

function addVisibleModifiers(name, modifiers) {
  if (modifiers.length === 0) return name;
  return `${name} (${modifiers.join(", ")})`;
}

function coreBaseCategory(baseName) {
  if (baseName === "Bay" || baseName === "Wildbay") return "BayOrWildbay";
  return baseName;
}

function classifyCrPrl(CrPrl) {
  // Cream/Pearl liegen auf demselben Locus: Cr / cr / pl
  // Cream: Crcr => 1, CrCr => 2, Crpl => 1 (wird später wie 2 benannt), sonst 0
  const cream = CrPrl === "CrCr" ? 2 : CrPrl === "Crcr" || CrPrl === "Crpl" ? 1 : 0;
  const pearl = CrPrl === "plpl" ? 2 : 0;
  const pearlCarrier = CrPrl === "crpl";
  const isCrpl = CrPrl === "Crpl";
  return { cream, pearl, pearlCarrier, isCrpl };
}

function computeColorName({ base, hasDun, hasChampagne, CrPrl }) {
  const baseCat = coreBaseCategory(base);
  const { cream, pearl, pearlCarrier, isCrpl } = classifyCrPrl(CrPrl);

  const chOn = hasChampagne;
  const dOn = hasDun;

  // Pearl (plpl)
  if (pearl === 2) {
    // Offene Kombi laut "~folgt~"
    if (dOn && chOn) return { name: "Unbekannt", pearlCarrier };

    if (dOn) {
      if (base === "Chestnut") return { name: "Apricot Dun", pearlCarrier };
      if (baseCat === "BayOrWildbay") return { name: "Pearl Bay Dun", pearlCarrier };
      if (base === "Black") return { name: "Pearl Black Dun", pearlCarrier };
      if (base === "Sealbrown") return { name: "Pearl Brown Dun", pearlCarrier };
    }
    if (chOn) {
      if (base === "Chestnut") return { name: "Gold Pearl", pearlCarrier };
      // ~folgt~: unknown
      return { name: "Unbekannt", pearlCarrier };
    }

    if (base === "Chestnut") return { name: "Apricot", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Pearl Bay", pearlCarrier };
    if (base === "Black") return { name: "Pearl Black", pearlCarrier };
    if (base === "Sealbrown") return { name: "Pearl Brown", pearlCarrier };
  }

  // Cream/Pearl Kombi (Crpl) wird wie CrCr benannt
  const creamLike = isCrpl ? 2 : cream;

  // D x Ch (ohne Cream/Pearl)
  if (dOn && chOn && creamLike === 0) {
    if (base === "Chestnut") return { name: "Gold Dun", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Amber Dun", pearlCarrier };
    if (base === "Black") return { name: "Champagne Dun", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sable Dun", pearlCarrier };
  }

  // D x Cr (Crcr)
  if (dOn && !chOn && creamLike === 1) {
    if (base === "Chestnut") return { name: "Dunalino", pearlCarrier };
    if (base === "Bay") return { name: "Dunskin", pearlCarrier };
    if (base === "Wildbay") return { name: "Wild Dunskin", pearlCarrier };
    if (base === "Black") return { name: "Smoky Grulla", pearlCarrier };
    if (base === "Sealbrown") return { name: "Smoky Brown Dun", pearlCarrier };
  }

  // D x CrCr / Crpl
  if (dOn && !chOn && creamLike === 2) {
    if (base === "Chestnut") return { name: "Cremello Dun", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Perlino Dun", pearlCarrier };
    if (base === "Black") return { name: "Smoky Cream Dun", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sealbrown Cream Dun", pearlCarrier };
  }

  // Ch x Cr (Crcr)
  if (!dOn && chOn && creamLike === 1) {
    if (base === "Chestnut") return { name: "Gold Cream", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Amber Cream", pearlCarrier };
    if (base === "Black") return { name: "Classic Cream", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sable Cream", pearlCarrier };
  }

  // Ch x CrCr / Crpl
  if (!dOn && chOn && creamLike === 2) {
    if (base === "Chestnut") return { name: "Cremello Champagne", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Perlino Champagne", pearlCarrier };
    if (base === "Black") return { name: "Smoky Cream Champagne", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sealbrown Cream Champagne", pearlCarrier };
  }

  // D x Ch x Cr (Crcr)
  if (dOn && chOn && creamLike === 1) {
    if (base === "Chestnut") return { name: "Gold Dun Cream", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Amber Dun Cream", pearlCarrier };
    if (base === "Black") return { name: "Classic Dun Cream", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sable Dun Cream", pearlCarrier };
  }

  // ~folgt~: D x Ch x CrCr / Crpl
  if (dOn && chOn && creamLike === 2) {
    return { name: "Unbekannt", pearlCarrier };
  }

  // Fallback: Single-Cream
  if (!dOn && !chOn && creamLike === 1) {
    if (base === "Chestnut") return { name: "Palomino", pearlCarrier };
    if (base === "Bay") return { name: "Buckskin", pearlCarrier };
    if (base === "Wildbay") return { name: "Wild Buckskin", pearlCarrier };
    if (base === "Black") return { name: "Smoky Black", pearlCarrier };
    if (base === "Sealbrown") return { name: "Smoky Brown", pearlCarrier };
  }

  // Fallback: Double-Cream / Crpl
  if (!dOn && !chOn && creamLike === 2) {
    if (base === "Chestnut") return { name: "Cremello", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Perlino", pearlCarrier };
    if (base === "Black") return { name: "Smoky Cream", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sealbrown Cream", pearlCarrier };
  }

  // Fallback: Champagne (ohne Dun/Cream)
  if (!dOn && chOn && creamLike === 0) {
    if (base === "Chestnut") return { name: "Gold Champagne", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Amber Champagne", pearlCarrier };
    if (base === "Black") return { name: "Classic Champagne", pearlCarrier };
    if (base === "Sealbrown") return { name: "Sable Champagne", pearlCarrier };
  }

  // Fallback: Dun (ohne Champagne/Cream)
  if (dOn && !chOn && creamLike === 0) {
    if (base === "Chestnut") return { name: "Red Dun", pearlCarrier };
    if (baseCat === "BayOrWildbay") return { name: "Classic Dun", pearlCarrier };
    if (base === "Black") return { name: "Grulla", pearlCarrier };
    if (base === "Sealbrown") return { name: "Brown Dun", pearlCarrier };
  }

  return { name: base, pearlCarrier };
}

function derivePhenotype({ E, A, G, CrPrl, Ch, Z, D, O, SPL, Fl, Sty, Ra }) {
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

  // 2) Farbnamen-Logik (D + Ch + Cr/pl)
  const hasDun = D !== "dd";
  const hasChampagne = Ch !== "chch";
  const { name: baseColorName, pearlCarrier } = computeColorName({ base, hasDun, hasChampagne, CrPrl });

  // 3) Silver
  const zCount = Z === "ZZ" ? 2 : Z === "Zz" ? 1 : 0;
  const isSilver = zCount > 0 && hasBlackPigment;
  let silverNote = null;
  const baseName = baseColorName;
  if (zCount > 0) {
    if (hasBlackPigment) {
      silverNote = zCount === 2 ? "Silver (homozygot)" : "Silver";
    } else {
      silverNote = "Silver (ohne Wirkung bei Fuchs)";
    }
  }

  // 4) Präfix/Suffix-Anhänge (Reihenfolge wie gewünscht)
  const prefixes = [];
  const suffixes = [];

  // Präfixe: Flaxen → Silver → Sooty
  if (Fl === "flfl" && base === "Chestnut" && CrPrl === "crcr") prefixes.push("Flaxen");
  if (isSilver) prefixes.push("Silver");
  if (Sty !== "stysty") prefixes.push("Sooty");

  // Suffixe: Rabicano → Overo/Splashed White/Pinto
  if (Ra !== "rara") suffixes.push("Rabicano");

  const hasOvero = O !== "oo";
  const hasSplash = SPL !== "splspl";
  if (hasOvero && hasSplash) suffixes.push("Pinto");
  else if (hasOvero) suffixes.push("Overo");
  else if (hasSplash) suffixes.push("Splashed White");

  const withAffixes = `${prefixes.join(" ")}${prefixes.length ? " " : ""}${baseName}${
    suffixes.length ? " " : ""
  }${suffixes.join(" ")}`.trim();

  // 8) Grey überschreibt (langfristig)
  const isGrey = G !== "gg";
  const tags = [];
  if (pearlCarrier) tags.push("Pearl Träger (crpl)");
  if (silverNote) tags.push(silverNote);
  if (isGrey) tags.push("Grey");

  // Grey: immer nur "Grey" + voller Name ohne Grey in Klammern
  const shown = isGrey ? `Grey (${withAffixes})` : withAffixes;
  const detail = null;
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

  // Warnungen (z.B. Letalität)
  const overoDist = distsByLocus.find((x) => x.locusKey === "O")?.dist;
  const lethalOO = overoDist ? (overoDist.get("OO") ?? 0) : 0;
  const warnings = [];
  if (lethalOO > 0) {
    warnings.push(
      `<div class="warn"><b>Warnung:</b> Erwartete Wahrscheinlichkeit für <b>OO</b> (Overo homozygot) ist ${roundPct(
        lethalOO
      )}% und gilt als letal.</div>`
    );
  }

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
    ${warnings.length ? `<div style="margin-top:10px">${warnings.join("")}</div>` : ""}
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
    O: "oo",
    SPL: "splspl",
    Fl: "FlFl",
    Sty: "stysty",
    Ra: "rara",
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

