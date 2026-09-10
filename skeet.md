# SKEET

Meiosis as a truth table. One diploid in, every gamete it can produce out,
laid out at once. This page covers what the skeet tab does, the genetics it
rests on, and how to hand a picked gamete to a physical DEF CON 34 badge.

Related:
- [`synthetic-population-genetics.md`](synthetic-population-genetics.md) - the genetics core (haploids, loci, meiosis, mating). Read that first for the math.
- [`gene-exchange.md`](gene-exchange.md) - the three-phase QR protocol skeet's Panel C rides on top of.
- [`eggstraction.md`](eggstraction.md) - the inverse map: many observed gametes from a physical badge, worked back to the diploid that produced them.
- [`plan-skeet.md`](plan-skeet.md) - design notes for the tab itself.
- [`nastea1/dc34-gamete/PROTOCOL.md`](https://github.com/nastea1/dc34-gamete/blob/main/PROTOCOL.md) - the QR wire format, byte by byte.

---

## 1. The Genome, In One Screen

A badge's genome is a **diploid**: an ordered pair of haploids, `(haplo0, haplo1)`.
Each **haploid** is nine one-byte loci, in this fixed order:

```
cd_period, cd_rate, cd_dir, sat, hue_ratedir, hue_base, hue_bound, chaser, nonlin
```

So a full genome is 2 x 9 = 18 bytes. That's the whole thing. See
[`synthetic-population-genetics.md`](synthetic-population-genetics.md) §1
for the per-locus role.

Skeet's Panel A shows that diploid in the same locus-row-major layout the
popsim / borg / vim gene inspectors use, with a `haplo0` column and a
`haplo1` column and no phenome column (skeet is about how gametes form, not
how they light up).

---

## 2. Mating, In One Paragraph

When two badges mate, each parent produces one **gamete** (a fresh 9-byte
haploid derived from that parent's diploid) and the child's diploid is
`(parent_A_gamete, parent_B_gamete)`. Mutation is applied to the gametes on
the way in. That's the whole mating process; the child badge type is
inherited from the receiver's slot 0. Full detail in
[`synthetic-population-genetics.md`](synthetic-population-genetics.md) §7.

Skeet does not mate. Skeet only enumerates one parent's gametes. Mating
happens on the physical badge, or in popsim / borg.

---

## 3. Haploid vs. Gamete

The words look interchangeable and technically the *type* is the same (nine
bytes), but they mean different things:

- A **haploid** is one of the two strands already sitting in a diploid.
  `haplo0` and `haplo1` are haploids. They are the badge's stored genome.
- A **gamete** is a new haploid *derived from a diploid* through meiosis:
  for each linkage group, coin-flip which parent haploid contributes.

A gamete is what you hand to a mate. A haploid is what you keep. The
special case where meiosis happens to pick every group from `haplo0` gives
you a gamete that is byte-identical to `haplo0`, but it is still the output
of a meiosis coin-flip, not a copy operation.

---

## 4. Linkage Groups: Why It's 32, Not 512

The nine loci do **not** segregate independently. They form five
**linkage groups**, each group sharing one coin flip in meiosis:

| group | bit | loci                              |
|-------|-----|-----------------------------------|
| A     | 0   | `cd_period, cd_rate, cd_dir`      |
| B     | 1   | `sat`                             |
| C     | 2   | `hue_ratedir, hue_base, hue_bound`|
| D     | 3   | `chaser`                          |
| E     | 4   | `nonlin`                          |

Within a group the loci co-segregate perfectly - they always come from the
same parent haploid in a given gamete. Across groups, the coin flips are
independent. This structure lives in `Diploid::meiosis()`; see
[`synthetic-population-genetics.md`](synthetic-population-genetics.md) §3
for the vectorized form and a warning about the easy mistake of drawing
nine coins instead of five.

Five coins => `2^5 = 32` distinct coin patterns. Every gamete a badge can
produce corresponds to exactly one of those 32 patterns. If the naive
"nine coins" model were right it would be `2^9 = 512`, and every locus
would look more diverse than it actually is. The same 32-not-512 fact is
what makes reconstructing a physical badge's diploid tractable from
observed gametes; see [`eggstraction.md`](eggstraction.md) §2.

A 5-bit pattern encodes one gamete: bit `g` = 0 means group `g` comes
from `haplo0`, bit = 1 means from `haplo1`. Skeet labels rows with the
5-bit binary string in MSB-left order (leftmost bit = group E =
`nonlin`, rightmost = group A = the `cd_*` block).

If the parent is homozygous at some group's loci (`haplo0` and `haplo1`
agree there), the two patterns that differ only in that group's bit
produce byte-identical gametes. Skeet shows both rows anyway. The view is
"all 32 mathematically possible patterns," not "unique byte outputs" -
the redundancy is data about how much of the parent's linkage structure
is homozygous.

---

## 5. What The Skeet Tab Does

**Panel A - parent diploid.** The badge you sent to skeet. Locus rows,
`haplo0` and `haplo1` columns, badge type in the header, `clear` button.
Populated by pressing `send to skeet` from any specimen inspector
(popsim, borg, vim gene). Sending is a snapshot; edits in vim gene
afterward do not retroactively update the skeet slot.

**Panel B - all 32 gametes.** 32 rows, one per coin pattern. Each row shows:

- the 5-bit pattern label,
- nine boxed locus cells filled with the resulting byte value, colored by
  source haploid so the linkage blocks are visible at a glance (cells
  belonging to the same linkage group always share a color in a given
  row), and
- a `+ Pick` button.

Row order is shuffled once when a badge is sent to skeet and held stable
while you compare, so the display doesn't jump under you between renders.
Sending a different badge triggers a fresh shuffle.

**Panel C - QR gene exchange.** Hidden until a gamete is picked. Then it
reveals a three-phase panel that mirrors vim gene's QR exchange 1:1,
minus the "gamete source" picker - here the source *is* the picked row.

Skeet does not run mutation, phenotype expression, or fertilization.
Mutation happens on the receiving badge when it accepts the sealed
gamete; that's what the third phase of the exchange is for. See
[`gene-exchange.md`](gene-exchange.md).

---

## 6. Using Skeet With A Physical Badge

You need: your DC34 badge (stock firmware is fine), this sim running in a
browser, and a phone or QR-reader app. No camera on the page.

1. **Load a diploid.** In popsim, borg, or vim gene, select or build the
   badge whose gametes you want to enumerate, and press `send to skeet`.
2. **Pick a gamete.** In Panel B, scan the 32 rows and press `+ Pick` on
   the one you want to hand to your physical badge. Panel C appears and
   loads that gamete as the seal payload. Picking a different row swaps
   the payload; picking the same row again unpicks and hides Panel C.
3. **Extract the nonce.** On your badge, display its phase-1 QR, scan it
   with your phone, and paste the decoded string into Panel C's phase-1
   text area. Press `extract nonce`. See
   [`gene-exchange.md`](gene-exchange.md) PHASE 1.
4. **Seal the gamete.** In Panel C phase 2, pick a `byte 15 :: badge type
   override` (defaults to the parent diploid's badge type; picking a type
   different from your physical badge's avoids the same-type inbreeding
   penalty) and press `seal`. A QR appears.
5. **Show your badge the QR.** Your badge scans it, treats it as a mate's
   gamete, mutates it on the way in, and updates its own diploid.
6. **(Optional) Round-trip check.** Panel C phase 3's `decode phase 2`
   button opens the sealed bytes back to a haploid preview under the
   same key and nonce, so you can confirm the recovered bytes match the
   gamete you picked.

The **k** row at the top of Panel C is the shared symmetric key `k0`
that authenticates the exchange. It's shared with vim gene: edit it on
either tab and both stay in sync. Default is the published `k0`; leave
it alone unless you're deliberately forking your own isolated
population. See [`gene-exchange.md`](gene-exchange.md) for the "k9000"
idea.

---

## 7. When To Use Skeet (vs. Vim Gene)

Vim gene edits a diploid and mints **one** gamete on demand - a random
draw from the parent's 32 possible gametes. Good when you don't care
which coin pattern you get, or when you want to send a specific strand
verbatim (`haplo0` or `haplo1`).

Skeet is read-only on the diploid and enumerates **all 32** at once.
Good when you want to:

- pick a specific coin pattern to hand to your badge (e.g. "give me the
  gamete where `chaser` came from `haplo1` but everything else came from
  `haplo0`"),
- see the parent's homozygosity structure directly (identical-looking
  rows across a bit-flip = that group is homozygous),
- inspect the full mating potential of a diploid you found interesting
  in popsim or borg without committing to any single draw.

The QR exchange payload is the same wire format either way; a physical
badge cannot tell whether a gamete was sampled or hand-picked.
