# ANCESTRY COMPOSITION FOR BADGES

## What Happens When You Point AncestryDNA's Patent At A DC34 Badge

Related: [synthetic-population-genetics.md](synthetic-population-genetics.md).

## Related repositories

- [`bunnie/dc34-api`](https://github.com/bunnie/dc34-api) - genetics core.
- [`bunnie/dc34-vault`](https://github.com/bunnie/dc34-vault) - exchange state machine.
- [`bunnie/dc34-console`](https://github.com/bunnie/dc34-console) - LED renderer.
- [`nastea1/dc34-gamete`](https://github.com/nastea1/dc34-gamete) - QR wire format.
- [U.S. Patent 12,626,778](https://patents.justia.com/patent/12626778) - "Accelerated hidden Markov models for genotype analysis" (issued May 12, 2026, priority Apr 13, 2022).

## The Rabbit Hole

While digging through DC34 genetics we realized the badge is doing, in
miniature, exactly the setup a real ancestry-inference pipeline works on:
diploid individuals with per-window linkage structure, labeled reference
distributions per subpopulation, and a mutation process with a known kernel.
So we went looking for how the real ones actually do it, found
[U.S. Patent 12,626,778](https://patents.justia.com/patent/12626778) - the
two-stage HMM AncestryDNA uses for ethnicity composition - and decided the
right thing to do was implement it on the toy.

The badge is genuinely a near-perfect target for this. Real ancestry
pipelines have to fake the founder distributions (they can't sample your
great-great-grandmother, so they use present-day people with well-documented
ancestry as a proxy - that's why their categories are big regions and not
villages). We don't have that problem: the founder distribution for each
badge type is literally written down in `Haploid::from_type`, one uniform
per locus per type. Ancestry inference on the badge starts from an
**analytically-known founder panel**, which is the modeling luxury a real
pipeline can only dream about.

This document is the "why and what" of the exercise, with the
implementation in §5.

---

## 0. What Even Is This Patent

Patent 12,626,778 is a **two-stage hidden Markov model** for inferring a
person's ancestry composition from a phased diploid genotype. The moving
parts:

- **Input**: a phased genotype - two haplotypes, one from each parent, both
  broken into windows of adjacent SNPs (claim 6).
- **Stage 1 (full-ethnicity HMM)**: an HMM whose hidden states are every
  supported ancestry label (~40 in the shipped AncestryDNA product) runs on
  each haplotype independently. Keep the labels that any window's decode
  actually visits, discard the rest (claims 1, 4).
- **Handoff**: combine the two per-haplotype candidate sets into a small
  candidate subset for this individual.
- **Stage 2 (simplified HMM)**: build a joint HMM whose hidden states are
  `(parent-1 label, parent-2 label, switch flag)` triples drawn only from
  the pruned candidate subset (claims 2, 3). The switch flag is a phasing
  swap between adjacent windows.
- **Decode**: run forward-backward on the joint HMM, count label
  frequencies along the posterior path, apply a min-threshold plus
  delta-rounding step so the displayed composition sums to exactly 100%
  (claims 7, 10). Display (claim 9).
- **Panel**: full-HMM transition probabilities are estimated from a
  labeled reference panel of individuals with known ancestry
  (claims 8, 12).

The whole reason for the two-stage structure is speed. Real ancestry
decodes have hundreds of labels and hundreds of thousands of SNP windows;
running the full label-cross-label joint at every window would be
prohibitive. Stage 1 is cheap because it's a 1D decode per haplotype;
stage 2 is expensive but only runs on the pruned lattice.

At DC34 scale - 8 badge types, 4 windows (see §2 for why 4 and not 5) -
the "acceleration" is irrelevant. Stage 2 finishes in microseconds either
way. We're implementing the two-stage architecture for **faithfulness to
the patent**, not for speed. This is a demonstration of what a real
ancestry pipeline is doing under the hood, run on a target small enough
to visualize every internal quantity.

---

## 1. Why The Badge Is A Near-Perfect Target

Three things a real ancestry pipeline has to fight, and the badge just...
hands them to us for free.

### 1.1 Labeled founder distributions in closed form

Real Ancestry can't sample your great-great-great-grandmother. What they
do is collect **modern-day people with well-documented recent ancestry**
(all four grandparents from region X) and use *those* as a proxy for the
founder distribution of region X. It works because regional allele
frequencies don't drift much across a few generations of low mixing, but
it's a real modeling assumption and it's why their labels are
"Scandinavian" rather than "your grandmother's hometown."

The badge tells us the founder distribution for every type analytically -
it's `Haploid::from_type` in `dc34-api/src/lib.rs`. Per-type uniform
priors over per-locus `u8` ranges (with the Goon `hue_base=0` and Uber
`hue_bound=255` spikes on top). The whole "panel" is one fixed
`8 × 9 × 256` tensor we can write down without simulating anything.

**One founder panel. Static. Analytic. Free.**

### 1.2 Free phasing

Real Ancestry has to run a global phasing algorithm to figure out which
allele at each heterozygous position sits on which chromosome (claim 11
of the patent gestures at this - SNP genotyping produces unordered pairs
and you have to reconstruct the ordering statistically). The badge
stores `haplo0` and `haplo1` explicitly at every step. We skip the
entire phasing subsystem.

### 1.3 Explicit linkage structure

The five linkage groups from meiosis (§3 of the pop-gen doc) are exactly
the windows the patent's claim 6 wants. Within a group, alleles
co-segregated on one meiosis coin flip. Across groups, the flips were
independent. The window boundaries **are** the recombination points, and
`P(switch)` between adjacent windows is exactly 1/2 - no biological
linkage-decay math to fake.

### 1.4 The extras that real ancestry doesn't get

- **The mutation kernel is fully specified.** The Gray-code XOR mutation
  process (§5.2 of the pop-gen doc) is a known 256×256 stochastic matrix
  `K_rate` per mutation rate level. The emission distribution at
  generation `t` is a closed-form convolution `π · K^t` - no need to fit
  a "how has this allele changed since the founder" model, we can just
  raise the matrix to the right power.
- **We can evolve the population forward.** No real Ancestry customer
  will ever watch their inferred ancestry drift over generations. We
  can, because the population sim is right there in the same
  `index.html`.

---

## 2. What Changes When You Downscale The Patent To A Badge

The patent's math applies directly, but three of its design choices are
sized for the human-scale problem and need to be re-thought when the
problem is a 9-locus toy.

### 2.1 Windows: we merge D and E

The natural windowing follows the five meiosis linkage groups:

```
W_A  = (cd_period, cd_rate, cd_dir)
W_B  = (sat)
W_C  = (hue_ratedir, hue_base, hue_bound)
W_D  = (chaser)
W_E  = (nonlin)
```

But the shipped `phenotype()` code (§4.1 of the pop-gen doc) has that
lovely asymmetric line:

```rust
nonlin: self.0[0].chaser.saturating_add(self.0[1].nonlin),
```

The displayed `phenotype.nonlin` couples `haplo0.chaser` with
`haplo1.nonlin` - one locus from group D on one haplotype, and one
locus from group E on the other haplotype. That's a cross-window,
cross-haplotype coupling, and it breaks the HMM's per-window emission
factorization if D and E are separate windows.

The fix is to **merge them**. Groups D and E become a single joint
window `W_DE` emitting `(chaser_0, chaser_1, nonlin_0, nonlin_1)` plus
the phenotype-derived `nonlin`. The delta constraint
`phenotype.nonlin = sat_add(haplo0.chaser, haplo1.nonlin)` sits inside
`W_DE`'s emission factor as a hard constraint. HMM structure preserved,
nonlin bug modeled exactly under the shipped code.

Final window set: `{W_A, W_B, W_C, W_DE}` - four windows, not five.

### 2.2 Transitions: no `ρ`, just `f_T(t)`

Real ancestry decodes have a nontrivial adjacent-window transition
matrix because SNP windows are physically close on a chromosome - the
same ancestry block often spans many windows before a recombination
event switches lineages. The patent's transition probabilities are
estimated from the reference panel to capture this.

At DC34 scale, adjacent windows are **entirely different linkage
groups**, each with its own independent meiosis coin. Under panmictic
mating, adjacent-window ancestry labels within a single haplotype are
independent after one generation. The transition matrix reduces to
rank 1:

```
P(T_{i+1} = T' | T_i = T) = f_{T'}(t)
```

where `f_T(t)` is the current population fraction of type `T`. That's
the entire transition structure. Simple, correct, and honest about the
fact that the panmictic-baseline mating policy has no linkage
autocorrelation to model.

Non-panmictic policies (heavy assortative mating, or heavy use of the
inbreeding-pass amplification) could produce nonzero adjacent-window
type correlations. If we ever ship a non-panmictic policy in the
evolve-forward loop and want the HMM to model it faithfully, we'll
measure the empirical correlation from the population buffer and plug
it back in. Until then, `f_T(t)` alone is the whole transition matrix,
and this is honest rather than window-dressing.

### 2.3 Emissions: convolution, not mixture

There is a tempting simplification: model the emission as a mixture
between the founder prior and the mutation-equilibrium distribution:

```
E_t(a) = (1 − μ(t)) · π(a) + μ(t) · U(a)
```

where `μ(t) = 1 − (1 − p_mutate)^t` and `U` is uniform. This is wrong
at any small `t`. It pretends every mutated allele jumps to uniform,
but the Gray-code kernel produces **near-neighbor** smearing - a
one-bit Gray flip is a small step in ordinary space, not a uniform
resample. At `t = 1` under Baseline, `E ≠ 0.75π + 0.25U`; it's the
convolution `π · K`, which stays concentrated near the founder range.

We commit to the correct model: `E_{T,ℓ,t} = π_{T,ℓ} · K_rate^t`, with
`K^t` precomputed by repeated squaring on a fixed `t`-grid. The full
emission tensor is `8 × 9 × 256 × 9 grid points ≈ 165k floats` per rate
level, sub-millisecond to look up. This distinction is a modeling
correctness issue, not a taste issue.

### 2.4 Hidden state semantics

The stage-1 hidden state at window `W` on haplotype `h` is the
**founder-generation badge type whose lineage the alleles in this
window descend from** through the meiosis-coin history. That's a
well-defined label per (haplotype, window) given the full pedigree,
and it matches the patent's ethnicity-label semantics exactly.

It is not "the type whose emission best matches the observation." The
HMM is a generative model over lineage, not a soft-classifier over
types. That distinction matters when interpreting the composition
report: "62% Human" means "62% of your window-haplotype slots trace
their lineage back to a founder-generation Human," not "62% of your
alleles look Human-ish."

### 2.5 The prior is a confound: admixture vs classification

The stage-1 posterior at a window factors into prior times likelihood:

```
P(T | alleles) ∝ f_T · E_T(alleles)
```

`f_T` is the population fraction; `E_T` is the emission. On a strongly
informative locus (hue, chaser), `E_T` spans orders of magnitude across
types, so the likelihood crushes the prior and `f_T` is irrelevant. On a
weak or overlapping locus (sat, cd_dir, cd_period, plus the two uniform
loci), `E_T` is nearly flat across types, so the posterior collapses to
`f_T`: every ambiguous window defaults to "Human 70%."

That default is a confound. It answers the wrong question. The same HMM
expresses two different questions:

1. **Classification** - "which population was this badge drawn from?"
   Here `f_T` is the correct prior: if 70% of badges are Human, a random
   badge really is more likely Human.

2. **Admixture** - "what fraction of this genome descends from each
   founder type?" Here `f_T` is noise. A badge that is genuinely 100%
   Uber founders should read ~100% Uber no matter how rare Uber badges
   are in the population.

Real ancestry tools face the same split. The reference panel is whatever
it is: dominated by the people who happened to send in DNA, mostly
European (the "Human" of this system), under-sampled for rarer
ancestries (Uber, CtfContest, Goon). That over-representation is not
made up and not wrong, but folding it into the prior biases every
uninformative segment toward the majority and flattens the interesting
minorities out of view.

The adjustment is to remove the `f_T` term for the admixture estimate:

```
P_admixture(T | a) ∝ E_T(a)        # f_T set uniform over the active types
```

Equivalently, subtract `log f_T` from the log-posterior, or run the
decode with a uniform prior. A power prior `f_T^γ` interpolates: `γ = 1`
is the population (classification) view, `γ = 0` is the unbiased
admixture view, and `0 < γ < 1` is shrinkage.

The badgecestry tab defaults to the admixture prior (`γ = 0`) and keeps
the population prior available as the classification view.

---

## 3. The Two-Stage Decode In Concrete Steps

Given a target diploid (18 bytes, `haplo0` and `haplo1`) and a
generation-`t` context:

**Stage 1 - full-type HMM, twice.**

Run forward-backward on the 8-state HMM once per haplotype, across the
4 windows. Emission at (type `T`, window `W`) is the product of
per-locus emissions `E_{T,ℓ,t}` over `ℓ ∈ W`. Transitions are rank-1
`f_T(t)`.

Compute per-window per-type marginal posteriors. For each haplotype,
sum posterior mass across windows per type, and **keep the top 3
types**. This is our claim-4 pruning rule: fixed `k=3`, deterministic,
no threshold to tune. Result: candidate sets `C₀` and `C₁`, each
exactly 3 types wide.

**Stage 2 - simplified joint HMM.**

Build the joint HMM whose hidden states are
`(T₀ ∈ C₀, T₁ ∈ C₁, switch ∈ {0,1})` triples. State space size:
`3 · 3 · 2 = 18` per window.

Emissions factor per locus per haplotype using the `switch` flag to
choose which parent-lineage emits which haplotype's allele. On `W_DE`,
add the phenotype delta factor:

```
δ(phenotype.nonlin = sat_add(a_{s,chaser}, a_{1-s,nonlin}))
```

which is a hard 0/1 constraint (the observed phenotype must be
consistent with the allele values at this state, under the shipped
asymmetric `phenotype()` code).

Transitions across windows use `P(switch) = 1/2` (each meiosis coin is
independent) and `f_T(t)` on the type labels.

Run forward-backward. Get per-window joint posteriors of shape
`(4, 3, 3, 2)`.

**Composition, and display.**

Marginalize the joint posterior over `switch` and over the *other*
haplotype's label to get per-window per-haplotype-per-type posteriors.
Sum across the 4 windows and both haplotypes. Normalize to a
composition vector. Apply the claim-10 rule: drop any type below a min
threshold (default 5%), redistribute the mass proportionally among
survivors, round to integer percentages summing to exactly 100.

---

## 4. The Glass Box: Watching The HMM Think

This is the part we care about most. Ancestry.com renders a pie chart
and a marketing story; the actual HMM is invisible to the user, and
that is *fine* for a company that just needs to sell the result. But
we're doing this because the pipeline is genuinely interesting, and
hiding the internals would defeat the whole reason for building the
thing.

`badgecestry` puts the composition pie at the top and then, immediately
below it, a permanently-rendered "How the model got here" panel with
**every intermediate quantity the decode produced**. Five sub-views,
top-to-bottom, in the order the algorithm produces them. Every one of
them recomputes live as you drag the generation slider.

### 4.1 Panel 1 - Emission distributions `E_{T,ℓ,t}(a)`

A grid of 256-wide row heatmaps, one per (type, locus) pair. Each cell
shows the emission distribution over `a ∈ 0..=255` at the current `t`.
Overlaid: a vertical tick per haplotype marking the user's observed
allele at that locus.

At `t = 0` these are narrow uniform boxes over the founder ranges - you
can literally see the §2 pop-gen-doc table rendered as pixels, with the
Goon `hue_base` and Uber `hue_bound` spikes visible. Drag `t` up and
the boxes smear, first at the edges (Gray-code near-neighbors leak
mass out one bit at a time) and then progressively wider. At large `t`
most rows collapse to a near-uniform gray - and *that*, visually, is
why ancestry inference loses discriminative power over generational
time.

Two rows always render identical across types regardless of `t`:
`cd_rate` and `hue_ratedir`. Their founder priors are uniform-on-full-
`u8` for every type. Convolving identical priors with the same kernel
yields identical outputs, so these loci contribute exactly zero
log-likelihood difference between types for all `t`. The visual
identical-ness is the explanation.

### 4.2 Panel 2 - Stage-1 per-window per-type posteriors

Two heatmap matrices, one per haplotype. Rows: the 4 windows
(`W_A, W_B, W_C, W_DE`). Columns: the 8 types. Cell color intensity =
posterior mass. Numeric value overlaid.

Reading a row across the 8 types shows which types are competing to
explain that window. A high-contrast row ("Uber dominates, everyone
else is at 0.05") means this window is informative. A near-uniform gray
row ("all types roughly equal") means this window can't tell types
apart - which happens naturally for the uniform-prior loci and happens
at high `t` for loci where the priors have smeared into each other.

This is what the patent's claim-4 language ("labels the haplotype's
decode actually visits") means in concrete pixels.

### 4.3 Panel 3 - Stage-1 → stage-2 handoff, showing `C₀` and `C₁`

A labeled block, prominently placed between the stage-1 matrices and
the stage-2 grid:

```
Stage 1 pruned to top-3:
  C₀ (haplo0 candidates):  Uber (0.42)   Human (0.31)   Other (0.14)
  C₁ (haplo1 candidates):  Human (0.38)  Uber (0.29)    Goon (0.18)

Pruned out:
  Community (0.06)   Village (0.04)   CtfContest (0.03)   None (0.01)
```

Every candidate has its summed-across-windows posterior mass chipped
next to it. Pruned types show up below the line with their masses too,
so you can see the margin by which each was cut. If the pruning is
close ("we barely kept Uber over Community"), that's visible; if it's
overwhelming, that's visible too.

This block is what makes the two-stage patent architecture *legible*.
The whole reason stage 2 has 18 states instead of 128 is right here on
screen: the model made a decision about which types to consider
jointly, and this is that decision.

### 4.4 Panel 4 - Stage-2 joint state posterior (the prominent one)

For each of the 4 windows, a `3 × 3 × 2` posterior grid, laid out as
two side-by-side 3×3 heatmaps (one for `switch=0`, one for
`switch=1`). Rows indexed by `C₀`, columns by `C₁`. Cell color =
posterior mass on that joint state at that window.

Three things this view lets you see directly:

- **Which parent-of-origin assignment the model prefers.** If the
  `switch=0` matrix has its mass concentrated and `switch=1` is
  nearly zero, phasing is unambiguous at this window. If both
  matrices look similar, phasing is ambiguous - meiosis really could
  have swapped the parental origin here and the data doesn't
  distinguish. The user gets to see that ambiguity, not just the
  marginalized-away answer.
- **Where the composition percentages come from.** Marginalize this
  grid across `switch` → per-haplotype per-type posterior. Sum across
  the 4 windows and both haplotypes → composition vector. The
  marginalization is shown below the grid as small chips/arrows so
  the user can follow the arithmetic from raw posterior to displayed
  percentage.
- **The nonlin bug's asymmetry, directly.** `W_DE`'s emission includes
  the phenotype delta factor from §2.1, which is a hard constraint -
  so `W_DE`'s joint posterior is typically sparser and more decisive
  than the other three. There's a "swap haplo0 ↔ haplo1" control
  button in the UI; hit it, and `W_DE`'s grid visibly shifts while
  the other three don't. That is the nonlin bug of §4.1 of the
  pop-gen doc, made visible as a direct observable rather than as a
  statistical claim.

### 4.5 Panel 5 - Forward/backward message trace

The under-the-hood view. For each stage-1 HMM (per haplotype), the α
and β message vectors at each window rendered as small 8-bar charts.
Arrayed as a 4-window × 2-message × 2-haplotype grid of small
multiples.

This is for the reader who actually knows HMMs and wants to audit the
forward-backward math. It goes at the bottom of the panel and casual
users can scroll past it - but it is *shown*, not hidden behind a
click. "Model internals" here means "on screen," not "toggleable."

### 4.6 Live coupling to the generation slider

All five panels recompute and re-render on every `input` event of the
generation-`t` slider. Cost budget per redraw: emission-tensor lookup
(`O(9 × 8)` 256-vector indexes) + two 4-window × 8-state forward-
backward passes + one 4-window × 18-state joint pass + rendering
~20k canvas cells. Sub-10ms in vanilla JS. If the visible bottleneck
turns out to be canvas draw, throttle with `requestAnimationFrame`
coalescing - but don't preemptively optimize.

Dragging the slider from `t=0` to `t=64` is, structurally, an
interactive demonstration of ancestry inference losing discriminative
power over generational time. Emission heatmaps smear, stage-1
posteriors flatten, the `C₀`/`C₁` candidate sets churn, stage-2
grids lose their concentration, and the composition pie rebalances
toward `f_T(t)`. The whole story - "generational time destroys
ancestry signal" - is a single mouse gesture.

### 4.7 Data attribution

Every widget in the internals panel has a `data-source` DOM attribute
naming which decoder module produced it
(`full_stage.forward_backward`, `simplified_stage.joint_posterior`,
etc.). If you inspect any bar in DevTools, you can trace it back to a
specific line in the JS port. This is how we operationalize "not a
black box" beyond just "we drew a chart" - every pixel is
attributable to a specific computation.

---

## 5. Implementation: Python reference, then JavaScript

badgecestry is the first feature here with a standalone reference
implementation, and the two copies are built in a deliberate order: the
algorithm is worked out in Python first, then ported to JavaScript for
the interactive tab.

### 5.1 Python (`hmm/`) - the source of truth

A small NumPy package that is the authoritative statement of the math.
Each module maps to one stage of §3, and the decoder stages return every
intermediate quantity the internals panel needs, not just the final
composition:

- `founder_panel.py` - the analytic `8 × 9 × 256` founder-prior tensor
  (§1.1), closed form, no simulation.
- `mutation_kernel.py` - the per-rate 256×256 Gray-code kernels and the
  `K^t` precompute on the fixed `t` grid (§2.3).
- `full_stage.py` - stage-1 forward-backward per haplotype; returns the
  α/β messages, per-window posteriors, and the top-k candidate set.
- `simplified_stage.py` - stage-2 joint decode over
  `C₀ × C₁ × {switch}`, with the `W_DE` phenotype delta factor baked in
  (§2.1).
- `compose.py` - the claim-10 min-threshold + delta-rounding display
  step.
- `popsim_bridge.py` - a test double for the JS popsim getter, plus
  `decode_composition()`, the end-to-end harness that chains the stages
  the way the JS `loadDiploid` path does.

That "return every intermediate" contract is the load-bearing idea: the
α/β messages, per-window posteriors, and joint posteriors all come back
as named fields, so each of the five glass-box panels in §4 consumes a
specific, unit-testable output instead of re-deriving anything.

### 5.2 JavaScript (`index.html`) - the UX port

The interactive `badgecestry` tab lives inline in the single-page app, no
separate file. Its functions mirror the Python modules one-to-one:
`buildFounderPanel`, `buildEmissionTensor`, `snapTToGrid`,
`decodeHaplotype`, `decodeDiploid`, `composePercentages`. On top of that
port sits the one layer with no Python twin - the canvas renderer for the
five panels (§4) and the `window.badgecestry` surface (`loadDiploid`,
`decodeComposition`) that the other tabs hand diploids into.

### 5.3 The order

Python first: get the algorithm right somewhere every module can be
unit-tested in isolation and the numbers inspected directly. JavaScript
second: port the validated functions into the browser for the live,
slider-driven UX. The NumPy package is the authority; when the two
disagree, the JS port is the bug.

---

## 6. What Falls Out For Free

Because the HMM is a static analyzer that runs on any diploid, and
because we already have a 10k-individual sim with policy-swappable
selection and mating, a fistful of experiments become one-line
changes:

1. **Ancestry-decay curves.** Sample a founder Uber diploid at `t=0`,
   evolve it under panmictic Baseline for `t ∈ {0, 10, 50, 200}`, plot
   the Uber posterior. Concrete question: how many generations until
   the composition-reported Uber fraction drops below the claim-10 min
   threshold? Prediction: `t ∈ [30, 100]` at Baseline. Run it and
   find out.
2. **Inbreeding-pass impact.** Same experiment with the inbreeding
   pass on vs. off. Does isolation preserve ancestry legibility or
   destroy it faster? Prediction: inbreeding *accelerates* ancestry
   loss for the isolated type, because the extra mutation pass smears
   the founder prior faster than same-type meiosis restores it.
3. **Nonlin bug likelihood ratio.** Run the same population twice
   under the shipped-bug and symmetric `phenotype()` variants,
   compute the composition marginal likelihood under each model, log
   the ratio. This is the quantitative version of §9 item 7 of the
   pop-gen doc: "the nonlin bug matters" becomes "the log-likelihood
   ratio at generation 50 under Baseline is X, summed across 10k
   diploids."
4. **Cross-type introgression detection.** At what mixing rate does
   the HMM start reliably calling an Uber-chaser allele in a Human
   background? Reliability threshold: the Uber composition fraction
   exceeds the min-threshold in ≥90% of 1000 Human-background
   diploids.
5. **Founder-panel sensitivity.** Perturb the §2 pop-gen-doc ranges
   by ±10% and rerun composition on a fixed test set. Types whose
   composition changes by >10 percentage points under this
   perturbation are structurally under-identified - a real result
   about the badge design, not just about the HMM.
6. **Kernel grid resolution.** Verify empirically that snapping `t`
   to the nearest grid point vs. linearly interpolating `K^t`
   between grid points changes the composition by <1 percentage point.
   If so, freeze the snap-to-nearest decision.

---

## 7. What This Isn't

It's worth being explicit about the boundaries.

- **This is not a replacement for a source audit.** The white paper's
  reading of the shipped code is what tells us the emission priors,
  the linkage structure, and the nonlin bug in the first place. The
  HMM sits on top of that audit, it doesn't replace it.
- **This is not a security tool.** The genetics crypto (AES-256-GCM-SIV
  with `k0`) is entirely separate from anything the HMM touches.
  Ancestry inference operates on the plaintext diploid post-decrypt;
  it does not participate in the exchange.
- **This is not a claim to know real-world author intent.** The
  patent's authors solved the human-scale problem. We're using their
  architecture on a toy where the acceleration doesn't matter and
  the founder distributions are analytic. The math applies; the
  motivations are ours.
- **This is not going to be accurate at generation 200 under
  Apocalyptic.** No ancestry pipeline is. When the founder priors
  have fully smeared into the mutation-equilibrium distribution,
  there's nothing left to infer. Panel 1 will show you that
  happening; the composition pie will collapse to `f_T(t)`. That's
  not a bug in the HMM, it's information-theoretic reality.

---

## 8. Provenance Note

The two-stage architecture, the window decomposition, the switch-flag
formulation, and the min-threshold delta-rounding display step are
transcribed from patent 12,626,778. The DC34-specific pieces (windows
as linkage groups, `P(switch) = 1/2` from meiosis coins, the
convolution emission model, the `W_DE` phenotype-factor treatment of
the nonlin bug, the founder-lineage state semantics) are our
extensions, forced by the shape of the badge target.

The badge source itself is transcribed from
`dc34-api/src/lib.rs`, `dc34-vault/`, and the derivations already in
[synthetic-population-genetics.md](synthetic-population-genetics.md).

Nothing in the patent constrains its architecture to human ethnicity.
The same math runs on any labeled multi-class prior over discrete
alleles with a linkage-group structure. DC34 badges are one such
target, and a delightful one for demonstrating what a real ancestry
pipeline is doing under the hood.
