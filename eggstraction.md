# EGGSTRACTION

Non-destructive genome readout of a physical DEF CON 34 badge, using nothing
but the stock gene-exchange protocol the badge already speaks. Point enough
minted nonces at the badge, collect enough of the sealed gametes it hands
back, and the underlying diploid falls out of the statistics.

This page is exclusively about physical badges. There is no simulated
counterpart. The eggstraction tab is a nonce factory, a paste box for the
badge's returned QRs, an evidence store for the observed gametes, and a
statistical reconstruction of the diploid that produced them. If you don't
have a real badge in your hand, the tab does nothing useful.

Related:
- [`synthetic-population-genetics.md`](synthetic-population-genetics.md) - haploids, loci, meiosis, linkage groups, mutation. The whole math core the tab inverts. Read §3 (meiosis / linkage groups) and §5 (mutation) first.
- [`skeet.md`](skeet.md) - the *forward* direction: one diploid, all 32 possible gametes enumerated at once. Eggstraction is the inverse.
- [`gene-exchange.md`](gene-exchange.md) - the two-QR protocol eggstraction rides on top of, from the perspective of the receiver.

---

## 1. What we are trying to do

You have a physical DC34 badge. You want to know its diploid - the two
haploids sitting inside it, `haplo0` and `haplo1`, nine bytes each. The
badge does not expose "print your genome." The only channel out is the
two-QR gene-exchange protocol, which reveals **one gamete per exchange**:
nine bytes, meiosis-sampled from the badge's diploid, then mutated.

The badge is otherwise unmodified. Stock firmware. No JTAG, no serial peek,
no debugger. The only asymmetry we get to exploit is that in the two-QR
protocol we (the sim) play the *receiver*: we generate the nonce, the badge
plays *responder* and hands us a sealed gamete. Under `k0` and our nonce, we
decrypt to nine bytes of that badge's genome, plus a badge-type tag in byte
15.

Repeat enough times and the badge is forced to reveal its full linkage
structure - not by any single exchange, but by the accumulated distribution
of what it emits. That's eggstraction.

Sister to skeet:

- **skeet** takes one diploid (in the sim) and enumerates every gamete it
  *could* produce - the forward map, `Diploid -> {32 gametes}`, all at once.
- **eggstraction** takes one physical badge (in the world) and enumerates
  every gamete it *does* produce across many exchanges - the inverse map,
  `{gametes observed} -> Diploid`, coverage-driven and mutation-aware.

---

## 2. Why 32, not 512: linkage groups

A DC34 haploid is nine one-byte loci, but they don't segregate
independently. They form **five linkage groups**
([synthetic §3](synthetic-population-genetics.md#3-meiosis-independent-assortment-but-with-linkage-groups)):

| group | bit | loci                              |
|-------|-----|-----------------------------------|
| A     | 0   | `cd_period, cd_rate, cd_dir`      |
| B     | 1   | `sat`                             |
| C     | 2   | `hue_ratedir, hue_base, hue_bound`|
| D     | 3   | `chaser`                          |
| E     | 4   | `nonlin`                          |

Meiosis is **five coin flips, not nine**. Each flip decides whether that
group's loci come from `haplo0` or `haplo1`, and every locus in the group
inherits from the same side. So every possible gamete a badge can ever
produce is one of `2^5 = 32` combinations, no more. A gamete is fully
identified by its five-bit *coin pattern* `bit4 bit3 bit2 bit1 bit0`.

A naive "each locus segregates independently" model would say `2^9 = 512`.
That is wrong for this genome, and the discrepancy is where the whole
reconstruction argument lives. Because groups A and C bundle three loci
under one coin each, the observed gametes are *not* independent nine-byte
draws - they are heavily correlated within each group. Any one gamete tells
us about all three `cd_*` loci at once, and about all three `hue_*` loci at
once.

That correlation is what makes reconstruction cheap: seeing one gamete is
seeing five group-slot values simultaneously. And within a group, the
co-occurrence structure of *which values appear together* lets us
reconstruct the two haploid tuples even after mutation has scrambled some
of the individual bytes. See §6 on group-tuple concordance.

**Note on what we cannot recover.** We can never distinguish "this is
`haplo0`" from "this is `haplo1`". Meiosis picks independently per group,
so the h0/h1 label is not observable through the protocol - swapping the
two per group produces the same distribution. What we recover is exactly
the equivalence class the badge stores: five unordered `{v0, v1}` pairs,
one per group. `(a, b)` and `(b, a)` are genotypically identical and
observationally indistinguishable.

---

## 3. The two-stage mutation model, and which stage we see

Every mating exchange in the DC34 system has **up to three** mutation
events. Eggstraction has to reason about exactly one of them and ignore the
other two. Getting this straight is the difference between a working
reconstruction and one that drifts by 30% every scan.

From [synthetic §5.4 and §6](synthetic-population-genetics.md#54-lock_rate-and-final_rate),
the three passes are:

1. **Sender-side sperm mutation (unconditional).** Before the badge seals
   its gamete into the phase-2 QR, it runs `mutate(gamete, final_rate)` on
   it. This happens *every* exchange, regardless of who the receiver is,
   whether the badge types match, or whether the sim is on the other end.
   **This is the mutation eggstraction sees.**
2. **Receiver-side egg mutation (unconditional, receiver-only).** The
   receiver mutates its own egg before joining sperm and egg into a child
   diploid. This doesn't touch the sperm on the wire, and in our case there
   is no receiver egg - we (the sim) are not producing a child. So this
   pass simply does not happen from our perspective.
3. **Inbreeding pass (conditional, same-type mating only).** If sperm and
   egg badge types match, the sperm gets an *extra* mutation pass on the
   receiver side at an elevated rate, and the egg is mutated at that
   elevated rate too. Again: receiver-side, doesn't touch what came off the
   wire. We are not accepting the sperm into any diploid, so this does not
   happen either.

**Consequence for eggstraction:** every gamete we observe has been through
**exactly one** mutation pass - the sender's, at the sender's `final_rate`
- and nothing else. There is no inbreeding penalty to model, no receiver
egg to think about, no compounding of rates. Just one clean mutation pass
on the way out of the badge, and then the seal.

**What rate is that pass at?** In principle any of `None`, `Baseline`,
`Elevated`, `Radioactive`, `Apocalyptic`. In practice, effectively always
Baseline. The badge's `mutation_param` is an integrator over physical
input - side-rocker reversals, accelerometer orientation changes - and it
decays back to Baseline whenever the badge sits still. To do an
eggstraction scan you hold the badge still, point its camera at the
laptop, wait for a QR, and scan a QR back into your phone. There is no way
to run that procedure without giving the meter plenty of time to decay.
Elevated / Radioactive / Apocalyptic require the badge to be actively
agitated during the entire exchange. In the wild, every gamete you scan
will be a Baseline draw.

The tab hardcodes this assumption:

```
EGG_MUT_ASSUMED_RATE = "Baseline"    // one Gray-bit flip per mutated locus, 25% per locus
```

If your target badge really is being shaken during scans, the reconstructor
will still converge, just more slowly and with more mass leaking to
gray-neighbor candidates. It is not currently exposed as a knob.

---

## 4. The mutation math itself: Gray code neighbors

The single most important structural fact about DC34 mutation is that it
happens in **Gray code space**, not binary. This changes everything about
what the "wrong" values around the true value look like, and it is the
lever that lets us do statistical reconstruction cheaply.

### 4.1 What actually happens to a mutated locus

From [synthetic §5.2](synthetic-population-genetics.md#52-how-big-a-mutation-is-or-why-gray-codes-are-ingenious-here):

```
if this_locus_mutates:
    new_allele = gray_decode( gray_encode(old_allele) ^ (bits << shift) )
```

where `bits` is a rate-dependent mask (`Baseline = 0x01`, i.e. one bit) and
`shift` is uniform in `{0..7}`, resampled per locus per pass.

Under Baseline, that means: encode the current byte value into its 8-bit
Gray code, flip exactly one bit at a uniformly random position, decode
back. The result is a **Gray-1 neighbor** of the original.

### 4.2 Why Gray codes matter for reconstruction

In *standard* binary, flipping bit `k` of a `u8` always changes the value
by exactly `±2^k`. A bit-3 flip is always `±8`, regardless of the current
value - uniform and boring.

In *Gray* code, adjacent codes differ by exactly one bit. So flipping one
bit in Gray-space is a *monotone-neighborhood* move in ordinary space: the
mutated value is typically nearby in the integer sense. Sometimes it is
`±1`, sometimes a longer nearby jump (e.g. Gray-flipping bit 7 of value
`127` lands on `128`, one step in ordinary space; Gray-flipping bit 0 of
value `128` lands on `129`, one step). Baseline mutations produce mostly
*small* value shifts.

The eight possible Gray-1 neighbors of a value `v` are its **mutation
cloud** under Baseline. Precomputed as LUTs (`EGG_GRAY_ENCODE`,
`EGG_GRAY_DECODE`), each locus's forward distribution under one Baseline
pass looks like:

```
P(observe = v)           = 1 - p_mutate            (0.75 under Baseline)
P(observe = gray_nb_i)   = p_mutate / 8            (0.03125 each, i in {0..7})
```

Total mass on the eight Gray-1 neighbors: `0.25 / 8 * 8 = 0.25`. Total mass
staying put: `0.75`.

### 4.3 The critical implication: values are not identities

Here is the mental shift that eggstraction requires. In skeet you can look
at a gamete row and read off `cd_period=3` as a *fact about the diploid*.
In eggstraction, seeing `3` at that locus in one scan is *evidence*, not
identity. It could be:

- the true value at one of the haploids (prob 0.75 per scan);
- a Gray-1 mutation of a neighbor (prob 0.25 / 8 per originating neighbor
  per scan);
- and after enough scans, any of these possibilities gets exposed by the
  ensemble frequency, not by any single row.

**You cannot take any observed byte at face value.** You have to treat the
per-locus scan history as a histogram over 256 bins and ask which
combination of two true values best explains that histogram.

### 4.4 Reading Gray-neighborhood as a distance metric

For reconstruction, define the *Gray-1 neighborhood* of a value `v` as the
set of nine values reachable in zero or one Baseline flips:

```
gray_nbrs_1(v) = { v } ∪ { gray_decode(gray_encode(v) ^ (1 << s)) : s ∈ {0..7} }
```

Every observed byte at a given locus lies in `gray_nbrs_1(true_value)`
with probability `1.0` under Baseline (each observation is either the true
value or one of its 8 Gray-1 neighbors). So the true value is **always in
the observed set union its own Gray-1 neighbors**, and enumerating
candidates for a locus is:

```
candidates(locus) = ⋃_{v : histogram[v] > 0} gray_nbrs_1(v)
```

This is the candidate enumeration the tab uses (`eggstrLocusCandidates`).
It is exhaustive under Baseline: no candidate outside this set can have
produced any of the observed values in a single mutation step.

### 4.5 Composing to a per-scan likelihood

Fold the meiosis coin in. Each scan draws from one of the badge's two true
values (`v0`, `v1`) with prob 1/2 each, then mutates:

```
P(observe = b | true = {v0, v1}) = 0.5 * P(b | v0) + 0.5 * P(b | v1)
```

where `P(b | v)` is the per-locus mutation likelihood - `0.75` if `b == v`,
`0.25/8` if `b` is a Gray-1 neighbor of `v`, `0` otherwise.

The tab precomputes `P(b | v)` as a 256x256 float table
(`EGG_MUT_LIKELIHOOD_BASELINE`), so scoring a candidate `{v0, v1}` pair
against a whole per-locus histogram is a fast tight loop over occupied
bins.

---

## 5. Confidence, not answers: how the per-locus estimator works

For each of the 9 loci independently:

1. Enumerate all candidate pairs `{v0, v1}` (v0 ≤ v1) over the
   Gray-1-extended candidate set for that locus.
2. For each pair, compute the log-likelihood of the observed histogram:

   ```
   LL(v0, v1) = Σ_b  H[b] * log( 0.5 * P(b|v0) + 0.5 * P(b|v1) )
   ```

3. Take the argmax as the point estimate: `bestPair = argmax LL(v0, v1)`.
4. Confidence is a **full softmax** posterior over the candidate pairs
   under a uniform prior:

   ```
   conf = exp(bestLL) / Σ_pair exp(LL_pair)
   ```

   or equivalently, in numerically stable form,

   ```
   conf = 1 / Σ_pair exp(LL_pair - bestLL)
   ```

This is a proper posterior probability. `conf = 0.95` really does mean
"under a uniform prior over candidate pairs, this pair is 95% of the
posterior mass." That is what we report to the user.

### 5.1 Why full softmax, not top-2 gap

An earlier version of this estimator computed the confidence as the softmax
of the top two candidates only. It behaved badly in one obvious case: when
the true pair is `{v, v}` (homozygous) and there is a Gray-1 neighbor `v'`
that has picked up a few mutated observations, the pairs

```
{v, v}          # homozygous truth
{v, v'}         # heterozygous, but v' is really just a mutation
```

both fit the histogram *very* well - `v` gets most of the mass, `v'` gets
a little. The top-2 gap between them is narrow, so the estimator would
report ~50% confidence in the homozygous call even though every *other*
alternative was strongly disfavored. That was misleading.

Full softmax fixes this. The homozygous pair still dominates the total
posterior mass because the other 200+ candidates are all far worse. A
close second-place gets damped by having only one close competitor among
many, not by being one of two.

The user-visible effect: confidence for "true hom vs adjacent-neighbor
het" stays high (~0.8–0.9) once you have enough scans, instead of pegging
at 0.5 forever.

### 5.2 What confidence you should demand

The tab defines "locked" as **per-locus confidence ≥ 0.85**, and a linkage
group as locked when all of its loci are locked. This is a display
threshold; it does not gate anything. Even below 0.85, the point estimate
is displayed - it is just tinted "provisional." Above 0.85, we show it as
"locked" and stop nagging.

For the coverage line at the top of the panel, the tab reports:

```
P(diploid known) ≈ ∏_{li=1..9} conf(li)
```

- a product across the nine loci. This is "and-across-loci" logic: the
whole diploid is only known if every locus is individually known. It is
conservative - locus confidences are not independent (linkage groups
couple them) - but the phasing step in §6 already exploits within-group
coupling separately, so treating the product here as a chained-AND lower
bound is reasonable.

Rule of thumb, in practice:

- **6 exchanges** ≈ 85% of groups locked on typical heterozygous diploids.
- **10 exchanges** ≈ 99% locked, and the coverage line usually says
  "diploid known" for a badge with no unusual structure.
- **15+ exchanges** required in the presence of a lot of homozygous
  linkage groups - see §7.

---

## 6. Making the groups agree: within-group phasing and tuple concordance

The per-locus estimator returns two values `{v0, v1}` for each of the 9
loci, but the `v0`/`v1` labels are per-locus arbitrary. Linkage group A has
three loci, each of which independently decides which of its two values
to call `v0` and which `v1`. There is no *a priori* reason those three
labelings agree with each other.

They *have* to agree, though. Under the badge's meiosis model, every
scan drew all three A-loci from the same parental haploid. So if a scan
matched `est.v0` at `cd_period`, it *should* also match `est.v0` at
`cd_rate` and `cd_dir` - modulo mutation. The tab does two things to
enforce this.

### 6.1 Phase alignment against a group anchor

For each multi-locus group (A and C), pick the group's first locus as the
anchor. For each other locus in the group, count how many scans have:

- same anchor-side and same locus-side (`aSide == bSide`),
- vs cross-sides (`aSide != bSide`).

If cross > same, the locus's labels are backwards - swap them. This is
`eggstrPhaseGroups()`. The output is: within each group, all loci agree on
what "parent 0" means and what "parent 1" means, so the reconstructed
`haplo0` row and `haplo1` row in the header display are internally
consistent per group.

### 6.2 Tuple concordance: reading the joint distribution

Per-locus estimation has a blind spot around Gray-neighbor pairs. Consider
group A with three loci, on a badge that is actually *homozygous* on A
with tuple `(3, 20, 1)`. Every clean scan lands `(3, 20, 1)`. But roughly
`1 - 0.75^3 = 58%` of scans will have at least one locus mutated, so the
observed 3-tuples include `(3, 20, 1)`, `(2, 20, 1)`, `(3, 21, 1)`,
`(3, 20, 0)`, and so on. Each per-locus histogram has a strong mode plus
some Gray-1 sprinkle.

The per-locus estimator, looking at `cd_period` alone, cannot definitively
distinguish "true value is `3`, mutation cloud is Gray-1 neighbors of `3`"
from "true values are `{3, 2}`, both real haploid values, and the badge
is heterozygous at this locus." Both explanations fit the marginal
histogram.

The way out is to look at the **joint distribution** over the whole group
tuple. Under mutation, `P(all 3 A-loci clean per scan) = 0.75^3 ≈ 0.42`.
So on a homozygous group, roughly 42% of scans should land the *exact*
tuple `(3, 20, 1)`, and the rest should scatter to a mutation cloud of
close-by 3-tuples that individually appear less often. On a heterozygous
group with tuples `(3, 20, 1)` and `(150, 90, 4)`, roughly 21% of scans
should land clean on each true tuple, and the rest scatter around both.

`eggstrGroupTupleConsensus()` implements this:

- Count exact tuples across all scans.
- If the top tuple has ≥ 35% of scans and the second tuple has ≤ 15%, call
  the group homozygous with the top tuple.
- If the top two tuples both have ≥ 20% of scans, call the group
  heterozygous with those two tuples.
- Otherwise, defer to the per-locus estimator.

This override is only applied to multi-locus groups (A, C) with at least 4
scans. It resolves the gray-neighbor blind spot cleanly, because the joint
distribution has 42% mass on true tuples that no single-locus mutation
cloud can reproduce.

For single-locus groups (B, D, E), the per-locus estimator is the whole
story - there is no joint to look at.

---

## 7. Homozygosity: the sticky case

A group is **homozygous** if `haplo0` and `haplo1` carry the same value
(or tuple) at that group. Homozygous groups are the hard case for
reconstruction because they always emit the same true value, so we never
get to see the "other allele" that would tell us we had a heterozygous
pair. The mutation cloud around a true homozygous value is *identical* to
what you would expect from a heterozygous pair `{v, v'}` where `v'` is a
Gray-1 neighbor of `v` picked up only occasionally.

### 7.1 The rut

If you look at, say, window C (the three hue loci) and the values keep
repeating the same tuple most of the time with some Gray-1 scatter, the
temptation is to write it off as "obviously homozygous, I have the value,
move on." Correct instinct, but the tab has to justify it with math or it
would also produce false-positive hom calls on true-het pairs where one
value happens to appear more often in a small sample.

The rut is: **the naive "two most frequent values" reconstruction fails
here.** Homozygous groups have exactly one dominant value; forcing a
second one gives you a Gray-1 neighbor of the true value, which is wrong.
The reconstructor has to be willing to *return the same value twice* -
`{v0 = v, v1 = v}` - as the best hypothesis, instead of forcing two
different values.

### 7.2 How the estimator handles it

Per-locus (§5), the candidate set includes both `{v, v}` (homozygous) and
every `{v, v'}` (heterozygous with a Gray neighbor). Each is scored by
its own log-likelihood against the histogram. Under enough scans, the
homozygous pair `{v, v}` outscores the near-hom heterozygous alternatives
because:

- `{v, v}` predicts 100% of the meiosis mass lands on `v` before mutation,
  giving `P(observe = v) = 0.75` and `P(observe = v')` = `0.25/8 = 0.03125`
  for each Gray neighbor.
- `{v, v'}` predicts 50% of mass on `v` and 50% on `v'`, giving something
  like `P(observe = v) ≈ 0.5 * 0.75 + 0.5 * 0.03125 = 0.39` and
  `P(observe = v') ≈ 0.39`.

If the observed histogram really has ~75% of mass on `v` and ~3% on each
of `v`'s Gray neighbors, `{v, v}` fits about 10x better per scan (in
log-likelihood terms) than any `{v, v'}`. After ~10 scans this becomes
decisive.

Under fewer scans, the estimator honestly returns lower confidence rather
than committing early. And in a group with multiple loci, the *tuple*
consensus of §6.2 gets there faster: even 4-6 scans concentrated on one
3-tuple with only mutation-shape scatter around it is a strong homozygous
signal.

### 7.3 Presentation

Homozygous per-locus estimates render with a **HOM tag** in the
reconstructed-diploid header row, and both `haplo0` and `haplo1` cells
show the same value in a "matches hom" swatch color. This is the display
rule: HOM means the estimator's argmax has `v0 == v1`.

Confidence is still reported. HOM with `conf = 0.6` means "the argmax is
homozygous but a heterozygous alternative is still plausible." HOM with
`conf = 0.95` means "we are confident this group is homozygous." Neither
locks - even 0.99 does not prove homozygosity, because the badge could in
principle cough up the missed second allele on exchange 47.

### 7.4 The all-homozygous badge

The pathological case: the badge is homozygous on all five groups
(`haplo0` and `haplo1` are byte-identical). Every gamete is a mutation
cloud around a single true haploid. The reconstruction estimator still
converges - all five per-locus argmaxes settle on their true values with
increasing confidence - but the "coupon-collector complete" story is
meaningless here. There are no second alleles to collect. The coverage line just reports the product of per-locus
confidences, and it climbs asymptotically toward 1.0 with more scans.

The tab does not pretend to certainty in this case. If you see all five
groups displaying HOM with high per-locus confidence, you have a
homozygous badge; if there is even one locus with a low-confidence HOM
call, more scans will either resolve it as high-confidence HOM or surface
a second allele.

---

## 8. Coverage estimator: putting it all together

The status line at the top of Panel B ("gamete accumulator") says
something like:

```
[ 8 exchanges :: 4/5 linkage groups locked :: P(diploid known) ≈ 0.71 ]
```

That is composed as follows:

- **`n exchanges`**: raw count of observed gametes in the log.
- **`K/5 linkage groups locked`**: count of groups whose loci all have
  per-locus confidence ≥ 0.85, after phasing and tuple-consensus overrides.
- **`P(diploid known)`**: product of per-locus confidences,
  `∏_{li=1..9} conf(li)`. Chained-AND lower bound; conservative.

You are "done" when the coverage line reads something like `5/5 groups
locked :: P(diploid known) ≈ 0.99`. In practice this happens at 10-15
exchanges for a normal heterozygous badge; more for homozygous groups.

The `send to vim gene` button at the top of the reconstructed-diploid
header is enabled once the reconstruction is coherent enough to hand off.
It materializes an *ordered* representative (`haplo0 = est.v0`,
`haplo1 = est.v1` per locus, with the anchor-phased ordering per group)
into vim gene's diploid editor for downstream inspection, phenotype
expression, or mating.

---

## 9. What we deliberately do not do

- **No mutation simulation.** Everything in Panel B is a raw meiotic draw
  from the physical badge, already mutated by the badge's sender-side
  pass, exactly as it came off the QR. The tab does not run the sim
  genetics core.
- **No k0 override to bypass authenticity.** If the badge is running a
  custom firmware with a different key (see [gene-exchange.md](gene-exchange.md)
  on `k9000`), we cannot open its seals. This is a hard failure surfaced as
  an error, not silent data.
- **No inbreeding penalty.** Doesn't apply - we are not a badge, we do not
  supply a "receiver type" for the badge to check, and even if we did, the
  inbreeding pass runs receiver-side on the sperm-plus-egg after we
  already have the sender-side gamete off the wire. It cannot affect our
  observations.
- **No fertilization / expression.** This tab does not phenotype the
  reconstructed diploid. Use `send to vim gene` for that; vim gene shows
  the LED pattern and the phenotype breakdown for any diploid you hand it.
- **No mixed-badge detection.** An earlier design flagged "impossible"
  observations (e.g. a scan whose values could not have come from the
  running two-allele hypothesis for some group) and warned about
  mid-scan badge switches. The mutation-aware estimator makes this
  ill-defined - novel values are *expected*, that is what mutation looks
  like - so the check is gone. If you switch badges mid-scan, the tab
  cheerfully mixes their observations and produces a nonsense
  reconstruction. Watch the label and use `reset` between different
  physical badges.

The tab is deliberately narrow: it is a physical-badge readout
appliance. Downstream analysis lives in the other tabs.

---

## 10. Using eggstraction with a physical badge

*Placeholder - to be filled in with screenshots and the click-by-click
procedure once the physical procedure is documented.*

The short version of the flow: mint a nonce in Panel C, show the QR to
your badge camera, let the badge produce and display its phase-2 sealed
QR on the OLED, scan that with your phone or a QR reader app, paste the
decoded base45 string into Panel C's phase-2 textarea, press
`open & record`. The gamete lands in the accumulator; the coverage line
updates. Batch: mint several nonces at once and paste multiple sealed
strings one per line into the phase-2 textarea.

Full walkthrough with badge photos, camera framing, and troubleshooting
for the common decode-failed cases will land here.

---

## 11. Why this is worth doing

Three concrete use cases.

1. **Extract YOUR badge's diploid, without a debugger, using only built-in
   functionality.** Without this tool, you can only obtain your diploid by
   dumping the firmware state (or by trying to infer it from the LED pattern,
   which, good luck with that. The badge only displays the phenotype, not the
   genotype.)
2. **Bring YOUR badge into the sim.** Once the diploid is reconstructed, one
   click sends it to vim gene, popsim, or borg. From there, your badge is a
   full first-class citizen of the sim! Mate it against hypothetical
   partners, run generations forward, compare its phenotype under
   different mutation rates, look for rare children. Eggstraction is the
   bridge from atoms to bits.
3. **Empirical check on the fairness of meiosis.** The badge's firmware
   coin is *supposed* to be fair, one flip per group, uniformly random.
   The coin-pattern coverage histogram at the bottom of Panel B is a
   direct visual audit: 32 bins, one per 5-bit coin pattern, counting how
   many scans hit that pattern. Over 50-100 exchanges you can literally
   see whether the RNG is fair. Consistent skew is evidence of a firmware
   issue worth reporting.

None of these require the sim to know anything about a specific badge in
advance. The badge just has to be willing to mate - which is what it does
anyway.
