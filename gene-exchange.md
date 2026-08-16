# GENE EXCHANGE

How two DEF CON 34 badges actually mate — the two-QR protocol, the crypto behind
it, the role of `k0`, and the "k9000" idea of forking your own genetically
isolated badge population. Ends with the recipe most people actually want:
build an ideal mate in `vim gene` and cross it into your own badge.

Related:
[`nastea1/dc34-gamete/PROTOCOL.md`](https://github.com/nastea1/dc34-gamete/blob/main/PROTOCOL.md) — the wire format spec, byte-by-byte.
[`bunnie/dc34-vault`](https://github.com/bunnie/dc34-vault) — the firmware side of the exchange.

---

## 1. The two-phase exchange

A DC34 mating is asymmetric: **only one badge's genome updates per exchange**.
Call that badge the **receiver** (its genome will change). Call the other the
**responder** (it gives up a gamete but is not itself modified).

There are exactly two QRs, in order:

```
      RECEIVER (B)                    RESPONDER (A)
      -----------                     -------------
  1.  generate nonce
      show QR ---------------------->  scan QR, extract nonce
                                       mint one haploid gamete via meiosis
                                       seal it under B's nonce + k0
                                  <---- show sealed QR
  2.  scan A's sealed QR
      open under own nonce + k0
      -> get A's haploid
      run own meiosis to get an egg
      combine: (A's haploid, B's egg) -> new diploid
      apply mutation (elevated if inbreeding)
      express phenotype -> new light pattern
      user chooses KEEP or REVERT
```

The scanner tells the phases apart only by testing whether the first 16 bytes
of the decoded payload equal `DC34_HEADER`. If they do, it's phase 1 (a nonce
being offered). Otherwise it's phase 2 (a sealed gamete arriving).

**Only B is modified.** A did not scan B's genome; A only scanned B's nonce.
For A to also get offspring, the roles swap and the whole thing runs again:
A publishes a nonce, B seals under it, A syngamies. Two full exchanges to
change both badges.

### What "syngamy" means here

In real biology, syngamy is sperm + egg → zygote. In this system the receiver
plays the role of egg-carrier:

- The **haploid A sent** becomes slot 0 of B's new diploid.
- **B's own fresh meiotic gamete** (drawn from B's current diploid at the
  moment of scan) becomes slot 1.

So the new B is `(A_haploid, B_own_gamete)` — one strand from each side. Then
per-locus mutation runs, phenotype is expressed via
`synthetic-population-genetics.md` §4's expression rules, and the ring updates.

### Byte 15 and the inbreeding pass

The sealed plaintext is 16 bytes:

```
plaintext[0..9]   = 9 gene bytes (the haploid)
plaintext[9..15]  = 0x00 * 6     (padding)
plaintext[15]     = sender's badge type
```

Byte 15 is load-bearing. If A's badge type equals B's own, B's firmware calls
it inbreeding and raises the mutation rate before applying it. See
`synthetic-population-genetics.md` §5 (`RATE_VALUE`, `RATE_BITS`). This is
what makes matings between different badge types produce cleaner offspring
than matings within a type.

The safe universal answer for byte 15 is `7` (None). No badge capable of
breeding has type None, so a `byte 15 = 7` gamete never trips the inbreeding
pass on any receiver.

---

## 2. The crypto, briefly

Only two primitives ship over the air:

- **base45** (RFC 9285) — how the raw bytes become text a QR code can carry
  in Alphanumeric mode. Alphanumeric mode is one QR version smaller than byte
  mode for the same payload, which means larger modules and a better chance
  of a scan across a phone camera.
- **AES-256-GCM-SIV** (RFC 8452) with `k0` as the key and B's phase-1 nonce
  as the 12-byte nonce. Only the phase-2 sealed gamete uses this; the phase-1
  nonce is transmitted in the clear (it's just a random challenge, header +
  12 bytes).

Both are implemented in JavaScript in the vim gene QR panel and in
[`gamete-workbench.html`](https://github.com/nastea1/dc34-gamete/blob/main/gamete-workbench.html);
neither depends on WebCrypto, so the whole thing works from a file opened
straight off disk.

**Why GCM-SIV specifically:** the "SIV" (Synthetic Initialization Vector)
variant is nonce-misuse resistant. Ordinary AES-GCM catastrophically fails
under nonce reuse — an attacker who sees two ciphertexts under the same
`(key, nonce)` pair can recover the authentication key and forge messages.
GCM-SIV downgrades gracefully instead: reusing a nonce merely leaks whether
the two plaintexts were identical. Given that the DC34 nonce is only 96 bits
(2⁻⁴⁸ birthday collision) and any given badge exchanges many gametes in a
weekend, SIV is the honest choice.

### Where `k0` sits

`k0` is a single global 32-byte key shared by **every** DC34 badge on Earth.
It is not per-device. It doesn't identify who you are. It just says "the
sender was some DC34 badge." Two consequences:

- One badge's copy unlocks the protocol for all of them. That is exactly the
  weakness [`nastea1/dc34badge`](https://github.com/nastea1/dc34badge)
  exploited to publish the value.
- There is no key rotation and no revocation. Authenticity in this scheme
  proves "some DC34 badge, somewhere," never *which* badge.

The badge firmware carries a hardcoded prefix of `sha256(k0)` — the eight
hex chars `dca9ea49` at `dc34-vault/src/main.rs:42` — used by each badge to
validate its own key material. The vim gene panel's "k" row runs the same
check when the input is the default value; a custom key is accepted with a
warning that the oracle didn't match.

The published `k0` value is:

```
7ad84ed0e00aec0499ede65615e1da517c0150230d2abc6ec7b566e621e740b3
```

It only protects the light gene exchange. The vault's stored credentials
sit under a separate per-device key and are not affected by `k0` leaking.

---

## 3. k9000 — forking your own genetically isolated population

`k0` is a parameter, not a constant. The whole exchange works under *any*
32-byte key; `k0` is just the one bunnie's firmware happens to accept. Pick
any other 32 bytes — call it `k9000` — and you have a fully working DC34-style
population that is **genetically incompatible with real DC34 badges**:

- A `k9000` seal cannot be opened under `k0`. GCM-SIV's tag will fail.
- A `k0` seal cannot be opened under `k9000`. Same reason.
- Everyone running your vim gene panel with `k9000` in the k slot can trade
  gametes with each other, run meiosis, express phenotypes — the whole
  simulation — and none of it touches or is touched by the DC34 population.

This is the useful bit: it gives you an isolated experimental colony without
having to reflash firmware. Handy for:

- Running a home-lab breeding experiment where you don't want your test
  mutants leaking into the conference population if you happen to point a
  camera at a real badge.
- Teaching the protocol at a workshop where everyone loads a common
  `k_workshop` and can mate freely without worrying about accidentally
  ingesting a stranger's genome.
- Modeling founder-effect / island-population dynamics where two colonies
  drift apart, and only occasionally does a "migrant" appear (someone who
  has both keys and manually shuttles a haploid across).

There is nothing to build. In the vim gene panel:

1. **k row → verify & save**: paste your own 64 hex chars. The oracle check
   will fail (that's expected — it's not `k0`), and the status line will say
   *custom k accepted; oracle not checked*.
2. Everyone in your isolated colony uses the same custom key. Store it
   anywhere you'd store a shared secret.
3. Trade phase-1 nonces and phase-2 sealed QRs exactly as before. Real DC34
   badges pointed at your QRs will produce `tag failure` and refuse the
   payload, and vice versa.

No physical hardware in a k9000 colony — this only makes sense if the
population lives inside vim gene panels (or something equivalent). Real
badges can't have their `k0` swapped without firmware surgery.

Note that `k9000` isn't cryptographically stronger than `k0`; it's just
scope-narrower. A leaked `k9000` compromises exactly your colony, in the
same way a leaked `k0` compromises the whole DC34 population.

---

## Addendum — Building an ideal mate and crossing it into your real badge

The rest of this document assumes two participants. But most people who found
this repo want a solo workflow: sculpt a haploid they like in `vim gene`, then
convince their own physical badge to accept it. Here is how, end to end.

The vim gene panel is the **responder** (A) in this scenario. Your physical
badge is the **receiver** (B). We're going to run one full exchange in which
the vim-gene-hosted "ideal mate" seals a gamete for your real badge.

Prerequisites: your real DC34 badge with a working camera, this sim served
somewhere with a screen big enough to hold a QR, and a way to point the
badge's camera at that screen.

### Step 0. Sculpt the mate you want (in vim gene)

Use the diploid editor rows to set haplo0 and haplo1 to whatever genes you
want the *mate* to have. Remember: your badge will get exactly one haploid
from this mate, and meiosis on the mate side will pick, per locus, from
haplo0 or haplo1 with a coin flip. If you want a deterministic gamete (no
coin flip), fill both strands identically — a homozygous mate always sends
the same haploid.

If you don't care about determinism, leave the strands different and use the
Row-3 `meiotic (seeded)` source with `use fixed meiotic seed` checked, so
you at least get reproducibility across runs.

### Step 1. Get your badge's nonce

The panel can't ask your badge for its nonce — the badge's camera and screen
are the only exchange interface it has. So:

1. On your physical badge, initiate a mating (press whatever button its
   firmware uses to start the exchange; the badge displays its phase-1 QR).
2. Scan that QR with a phone camera or QR reader app. You'll get a base45
   string starting with the characters that decode to `DC34_HEADER`.
3. Paste that base45 into the vim gene panel's **decode** row (Row 4) and
   press *decode*. The status line will read *phase 1 :: nonce = XX XX ...*
   with 24 hex chars.
4. Copy those 24 hex chars.

If you can't easily paste the base45 into decode — say you scanned with a
phone that's not on the same clipboard — the phone's QR reader will show the
raw string. Type it in, or use any online base45 decoder to get the raw
bytes, and skip straight to the hex nonce (bytes 16..28 of the decoded
payload).

### Step 2. Feed the nonce into Phase 2 and seal

Back in the vim gene panel:

1. Open the **phase 2 :: sealed gamete** row.
2. Paste the 24 hex chars into *their nonce*. (If you decoded via the panel's
   Row 4, this field auto-fills — but only when you use Row 4's own nonce,
   which is the phase-1 QR case, not this scenario. For a real badge, paste
   by hand.)
3. Pick a gamete source. If you built a homozygous mate in Step 0, all three
   (`haplo0`, `haplo1`, `meiotic`) produce the same 9 bytes and it doesn't
   matter. Otherwise, pick `meiotic (seeded)` for reproducibility or
   `meiotic` for a fresh coin-flip gamete each seal.
4. **byte 15 :: badge type override**: pick a badge type that is **not** the
   same as your physical badge's own type. If your badge is Human, pick
   anything except Human. The safe universal choice is `7 :: None`, which
   never trips inbreeding on any receiver — but if your goal is a *clean*
   cross (Baseline mutation rate, one bit-flip per gene on average), picking
   any non-matching real type works and is more "in the world."
   Picking your own type on purpose is how you'd deliberately trigger the
   elevated (~39%) inbreeding mutation pass, if you wanted to see that
   effect on your badge.
5. Press **seal**. The QR renders under the nonce you pasted, with the
   plaintext bytes shown so you can double-check what's about to go over.

### Step 3. Show the QR to your badge's camera

Point your physical badge's camera at the QR you just generated in Row 3.
The badge is still in phase-2-waiting state from Step 1. It scans, opens
under its own nonce + `k0`, gets your ideal mate's haploid, runs its own
meiosis, syngamies, mutates, and expresses. Its ring updates.

If your badge shows a KEEP / REVERT prompt: this is where you decide whether
the mating "took." Keep it if the phenotype is what you wanted; revert if
the receiver-side meiosis diluted things too much or the mutation pass
scrambled a locus you cared about.

### Realistic caveats

- **You control half the resulting diploid, not the whole thing.** Your
  badge's own diploid still contributes one haploid via its own meiosis.
  A wide hue band from the mate lands in a single exchange (both `hue_base`
  min and `hue_bound` max rules pull toward the extremes), but a rare
  `chaser < 88` value from the mate usually does not, because your badge's
  own egg almost certainly has `chaser >= 88` and the `chaser`-expression
  rule is a saturating add. Expect several exchanges to fix a rare trait.
- **Mutation can undo your careful sculpting.** Even Baseline (~25% per
  locus) is a lot. If you want a specific gene to survive, prefer sending
  it as a homozygous haploid *and* pick a byte-15 that avoids the elevated
  inbreeding pass.
- **The panel doesn't emulate the receiver.** Row 4's *decode* shows the
  haploid your mate would send, not the diploid your badge would end up
  with after syngamy + mutation. To model that side too, use the popsim or
  borg tab, seed both parents into a two-individual population, and run one
  generation.

### If you want to do it entirely in-browser (no physical badge)

Same panel, roles played by clicks:

1. Row 2: generate a nonce for "receiver B" (this is just the panel
   pretending to be B).
2. Row 3: paste that same nonce into *their nonce* (auto-filled), pick your
   mate's source, seal.
3. Row 4: paste the Row-3 base45 output into *paste base45*, and paste the
   Row-2 nonce hex into *nonce*, then *decode*. The animated preview under
   the info block is the haploid A would have sent — homozygously expressed
   so you can see what it "looks like" if it fully dominated. It is not
   the syngamy result; for the real syngamy phenotype you'd combine that
   haploid with a haploid from your own diploid and run `phenotype()` over
   the pair, which the popsim tab already does at scale.
