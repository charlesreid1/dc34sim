# FAQ

- [What is this?](#what-is-this)
- [Using the simulator](#using-the-simulator)
- [Badgecestry](#badgecestry)
- [Your physical badge](#your-physical-badge)
- [Under the hood](#under-the-hood)
- [Related repositories](#related-repositories)

---

## What is this?

A browser simulation of the DEF CON 34 badge's population-genetics system:
10,000+ virtual badges evolving in your tab, plus a widget that can talk to
your physical DC34 badge over the QR exchange. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md).

It's a CRT-styled terminal with four tabs across the top. Each is a different
lens on the same genetics core:

- **popsim** - evolve a whole population of badges. Click any cell to inspect it.
- **vim gene** - edit one badge's genome (a diploid: 2 x 9 bytes) and mint QRs a real badge will accept.
- **borg** - same as popsim, but every mating uses one frozen partner, ROOT.
- **badgecestry** - decode one badge's ancestry composition with a two-stage hidden Markov model.

<!-- SCREENSHOT: full-page view with the four tabs circled -->

### Do I need a camera?

No. Nothing on the page needs a camera, and everything is click-driven. To
bring in a physical badge, you read the badge's screen with your own phone or
a QR reader app, then paste the text into the widget. See
[Your physical badge](#your-physical-badge) below.

## Using the simulator

### What's the difference between popsim and borg?

popsim simulates a free-mixing population of badges. Each generation, any
badge mates with any other at random.

borg runs the same simulator as one big colony: every badge mates with the
same fixed partner, ROOT, a frozen diploid. It's the social-insect colony
setup, one individual as the mate for everyone, so ROOT's genome spreads
through the colony until, without mutation, every badge is ROOT. See
[`borg-genetics.md`](borg-genetics.md).

### How do I evolve a population?

1. Open the **popsim** tab.
2. In the **control** panel (right side), set a population size and mutation rate.
3. Press **run** to evolve continuously, or **step +1** / **+10** / **+100** to advance one generation at a time.

See [`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 5 (mutation) and section 7 (mating).

### How do I select and inspect a simulated badge?

popsim and borg simulate badges (and the firmware running them) entirely in
the browser, so you can see and control the exact state of any individual.
On either tab, click a cell in the population grid. Its genome and phenotype
show up in the specimen inspector right below the grid.

The inspector is the hub between tabs: any badge you select can be sent
straight to **vim gene** to edit, or to **badgecestry** to decode. So the
whole pipeline is: pull a gene out of a popsim run, send it to vim gene, edit
it, and borg draws its frozen ROOT from the vim gene state.

### What do the allele-frequency histograms mean?

On the **popsim** tab (or **borg**), the **allele frequency spectra** panel
sits below the specimen inspector. Each histogram there is the live
distribution of one phenotype - hue, saturation, chaser, and so on - across
the whole population at this moment. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 9, "What To Measure".

### How do I send a specimen to vim gene?

The specimen has to be selected first, in either popsim or borg:

1. Open the **popsim** tab (or **borg**).
2. Click a cell in the population grid to select a badge. Its details appear in the specimen inspector below the grid.
3. In the inspector, press **send to vim gene**.
4. Switch to the **vim gene** tab. The badge is loaded into the diploid editor, ready to edit or export as a QR.

See [`gene-exchange.md`](gene-exchange.md) PHASE 0.

### How do I swap haplo0 and haplo1, and why would I?

In the **vim gene** editor, the **swap haplo0 ↔ haplo1** button exchanges the
two strands. It exists to expose the shipped `nonlin` asymmetry: the same two
haploids in opposite roles usually produce a different displayed phenotype.
See [`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 4.1.

## Badgecestry

### What is badgecestry?

It decodes a single badge's ancestry composition - which founder badge types
its genome descends from - using the two-stage hidden Markov model from U.S.
Patent 12,626,778 (the technique AncestryDNA uses). See
[`badgecestry.md`](badgecestry.md).

### How do I load my own gamete? (jizz in a cup)

In the **badgecestry** tab, use the **load from gamete** panel (top left):

1. Press **generate nonce**, then scan that QR with your badge.
2. The badge seals its gamete under that nonce and shows a second QR.
3. Scan the second QR and paste its text into the box.

Leave the box empty to use the provided genome.

This is where a physical badge differs from the simulated ones. popsim and
borg simulate badges in the browser, so you can inspect a whole diploid - both
haploids - directly. A real badge has no way to hand over its full genome.
There's no "here are my two strands" command. It can only mint a gamete: a
random new haploid assembled from those two strands. So loading your own
gamete is asking the badge to jizz in the cup and hand back one randomized
strand, not dump its whole genome.

<!-- SCREENSHOT: the load-from-gamete panel with the two QRs -->

### What does the "generations since founding" slider do?

It sets t, the number of generations of mutation since the founding
population. Drag it up and the ancestry signal decays. See
[`badgecestry.md`](badgecestry.md) section 4.6.

## Your physical badge

### Do I need my badge to use the widget?

No. Only the QR exchange needs it. Simulation, editing, and badgecestry all
run in the browser with nothing but the page.

### What is k0, and where do I get it?

k0 is the published default badge secret that encrypts the QR exchange. The
widget ships with it pre-loaded. See [`gene-exchange.md`](gene-exchange.md)
"Where k0 sits".

### How do I generate a nonce?

On the receiving badge (the one that will evolve):

1. On your badge, press left or right until it shows its phase-1 QR.
2. Scan that QR with your phone and copy the text it decodes to.
3. In the **vim gene** tab, open the **QR gene exchange** panel, phase 1, and paste the string.
4. Press **extract nonce**.

See [`gene-exchange.md`](gene-exchange.md) PHASE 1.

### How do I seal a gamete?

In the **QR gene exchange** panel, phase 2: choose a gamete source and press
**seal**. It produces a QR your badge will accept as a mate. See
[`gene-exchange.md`](gene-exchange.md) PHASE 2.

### What do I scan, and with what?

Your badge scans the QRs the widget draws. You scan your badge's screen with
your own phone and paste the text into the widget. There is no camera on the
page. See [`gene-exchange.md`](gene-exchange.md) PHASE 3.

<!-- SCREENSHOT: phone scanning the widget QR / badge screen -->

## Under the hood

Quick definitions only. Each term links to the document with the real math.

### What is a haploid / haplo0 / haplo1?

A haploid is one strand of the genome (9 bytes). A diploid is two strands,
haplo0 and haplo1, and every badge is diploid. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 1.

### What is the genome / "the nine loci"?

Nine one-byte loci (each 0-255) make up the genome. The phenotype you see on
the LEDs is computed from the two haploids. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md)
sections 1 and 4.

### What is mutation?

A per-locus chance to flip bits in a gamete (Gray-coded), set by the mutation
meter from Baseline up to Apocalyptic. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 5.

### What is the "inbreeding pass"?

When two mating badges share a badge type, both gametes take an extra
mutation pass - an assortative-mating penalty. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 6.

### What is the "nonlin" anomaly?

The shipped firmware expresses the nonlin phenotype asymmetrically (one
parent's chaser allele shows up in the other's nonlin). It shouldn't work, but
it ships anyway. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md)
section 4.1.

### What does "ROOT" mean?

ROOT is a frozen diploid used as the fixed partner in every borg mating, so
its genome steadily spreads through the colony. See
[`borg-genetics.md`](borg-genetics.md).

## Related repositories

- [`bunnie/dc34-api`](https://github.com/bunnie/dc34-api) - genetics core: `Haploid`, `Diploid`, `BadgeType`, `phenotype`, `meiosis`, `mutate`.
- [`bunnie/dc34-vault`](https://github.com/bunnie/dc34-vault) - exchange state machine and the k0 oracle.
- [`bunnie/dc34-console`](https://github.com/bunnie/dc34-console) - LED renderer.
- [`nastea1/dc34-gamete`](https://github.com/nastea1/dc34-gamete) - QR wire format (`PROTOCOL.md`), published k0, and the reference workbench.
- [U.S. Patent 12,626,778](https://patents.justia.com/patent/12626778) - "Accelerated hidden Markov models for genotype analysis".
