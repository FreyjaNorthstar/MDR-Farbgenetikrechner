const LOCI = [
  {
    key: "E",
    label: "Extension (E/e)",
    alleleOrder: ["E", "e"],
    genotypes: ["EE", "Ee", "ee"],
  },
  {
    key: "A",
    label: "Agouti (A/a)",
    alleleOrder: ["A", "a"],
    genotypes: ["AA", "Aa", "aa"],
  },
  {
    key: "G",
    label: "Grey (G/g)",
    alleleOrder: ["G", "g"],
    genotypes: ["GG", "Gg", "gg"],
  },
  {
    key: "Cr",
    label: "Cream (Cr/cr)",
    alleleOrder: ["Cr", "cr"],
    genotypes: ["CrCr", "Crcr", "crcr"],
  },
  {
    key: "D",
    label: "Dun (D/d)",
    alleleOrder: ["D", "d"],
    genotypes: ["DD", "Dd", "dd"],
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

function derivePhenotype({ E, A, G, Cr, D }) {
  // Sehr vereinfachte Ableitung, aber konsistent und erweiterbar.
  // 1) Basisfarbe via E/A
  const hasBlackPigment = E !== "ee";
  let base;
  if (!hasBlackPigment) base = "Fuchs";
  else base = A === "aa" ? "Rappe" : "Brauner";

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

  // 3) Dun
  const isDun = D !== "dd";
  let dunNote = null;
  let withDun = baseCream;
  if (isDun) {
    dunNote = "Dun";
    // grobe deutsche Bezeichnung
    if (baseCream === "Fuchs") withDun = "Fuchsfalbe (Red Dun)";
    else if (baseCream === "Brauner") withDun = "Falbe (Bay Dun)";
    else if (baseCream === "Rappe") withDun = "Mausfalbe (Grullo)";
    else withDun = `${baseCream} (Dun)`;
  }

  // 4) Grey überschreibt (langfristig)
  const isGrey = G !== "gg";
  const tags = [];
  if (creamNote) tags.push(creamNote);
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
    Cr: "crcr",
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

