"""The templating gate.

Four different pinned template chains sit in one soup of free monomers (equal numbers of
each letter, low density). Pump on. New backbone bonds form only by the type-blind,
held-gated ligation rule. Every ligation is logged, and for each new monomer we note the
template letter it sits against (if any template monomer is within 2.8 cells). Score:
how often the new monomer is the complement of that template letter. Chance is about 25%.

Round 2: the longest chain that formed at a template in round 1 becomes the template in a
fresh soup (four copies of it). Its copies should be ITS complement, which spells the
original again.

The reference for the score is chance: if linking at a template ignored letters, a new
monomer would be the complement about 25% of the time (equal numbers of each letter).

Arms (same seeds):
  receiver   the pumped-ear world (A+T = C+G = 2*Omega)
  null       an alphabet in which no two letters sum to within 0.04 of 2*Omega
             (A .02, C .03, G .04, T .05). NOT a clean null: at these low frequencies the
             phase sensitivity (~1/w) is large enough that each monomer's own sidebands
             give a second-order same-letter coupling. Kept for the record only.
  medium     the pump drives the field's mass term instead of the monomers' ears
  pump_off   pump switched off

    python experiments/gate.py --seeds 3 --round2
"""
import argparse, json, os, sys, time
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from chains import Config, World
from chains.analysis import complement_map, nearest, chain_report, duplex_report

TEMPLATES = ["GATTACAG", "CTGACCTA", "TCAGGATC", "AGCTTCGA"]
NULL_ALPHABET = {"A": 0.02, "C": 0.03, "G": 0.04, "T": 0.05}


def arm_config(arm, seed, N):
    cfg = Config(N=N, seed=seed)
    if arm == "pump_off":
        cfg.pump_on = False
    elif arm == "null":
        cfg.alphabet = dict(NULL_ALPHABET)
    elif arm == "medium":
        cfg.mode = "medium"
    return cfg


def run_round(cfg, templates, n_each, steps, log_every=5000, label=""):
    w = World(cfg)
    N = cfg.N
    centres = [(N / 4, N / 4), (3 * N / 4, N / 4), (N / 4, 3 * N / 4), (3 * N / 4, 3 * N / 4)]
    tids = []
    for seq, (cx, cy) in zip(templates, centres):
        ids = w.add_chain(seq, (cx - (len(seq) - 1) * cfg.d0 / 2, cy))
        tids += ids
    w.pinned[tids] = True
    free = np.repeat(np.arange(len(w.letters)), n_each)
    w.rng.shuffle(free)
    w.add_monomers(free)
    scoring = complement_map(Config())          # A-T, C-G, by the default physics
    seen = 0
    events = []          # (new letter, template letter it sits against)
    at_template = set()
    t0 = time.time()
    for t in range(steps):
        w.step()
        if len(w.events) > seen:
            new = w.events[seen:]
            seen = len(w.events)
            for (_, i, j) in new:
                part, _ = nearest(w, [i, j], tids, 2.8)
                for m, p in zip((i, j), part):
                    if p < 0:
                        continue
                    events.append((w.letters[int(w.typ[m])], w.letters[int(w.typ[p])]))
                    at_template.add(int(m))
        if log_every and (t + 1) % log_every == 0:
            print(f"  {label} t={t+1} ligations={w.ligations} template-site events={len(events)} "
                  f"({time.time()-t0:.0f}s)", flush=True)
    made = []
    for path in w.chains(min_len=2):
        if set(path) & set(tids):
            continue
        if set(path) & at_template:
            made.append(w.sequence(path))
    soup = chain_report(w, exclude=tids)
    dup = duplex_report(w, scoring, exclude=tids)
    ok = sum(scoring.get(b) == a for a, b in events)
    return {"ligations": w.ligations, "site_events": len(events), "site_correct": ok,
            "site_fidelity": ok / len(events) if events else None,
            "pair_counts": dict(Counter(f"{b}->{a}" for a, b in events)),
            "made_at_template": sorted(made, key=len, reverse=True),
            "soup_chains": len(soup), "soup_len_max": max([c["len"] for c in soup], default=0),
            "duplex": dup}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--steps", type=int, default=40000)
    ap.add_argument("--N", type=int, default=128)
    ap.add_argument("--n_each", type=int, default=53)     # same density as 30 per letter on 96x96
    ap.add_argument("--arms", nargs="+", default=["receiver", "pump_off", "medium"])
    ap.add_argument("--round2", action="store_true")
    ap.add_argument("--out", default="results/gate.json")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    res = {}
    for arm in a.arms:
        res[arm] = []
        for s in range(a.seeds):
            cfg = arm_config(arm, s, a.N)
            r = run_round(cfg, TEMPLATES, a.n_each, a.steps, label=f"{arm}/s{s}")
            r["seed"] = s
            if a.round2 and arm == "receiver" and r["made_at_template"]:
                seq2 = r["made_at_template"][0]
                cfg2 = arm_config(arm, 100 + s, a.N)
                r2 = run_round(cfg2, [seq2] * 4, a.n_each, a.steps, label=f"round2/s{s}")
                r2["template"] = seq2
                r["round2"] = r2
            res[arm].append(r)
            print(f"{arm} seed {s}: ligations {r['ligations']}, template events {r['site_events']}, "
                  f"complementary {r['site_correct']} ({r['site_fidelity']}), made {r['made_at_template'][:6]}, "
                  f"soup longest {r['soup_len_max']}" +
                  (f" | round2 on {r['round2']['template']}: {r['round2']['site_correct']}/{r['round2']['site_events']}"
                   if "round2" in r else ""), flush=True)
            json.dump(res, open(a.out, "w"), indent=1)
    print("\nSUMMARY (chance ~0.25)")
    for arm, rs in res.items():
        n = sum(r["site_events"] for r in rs)
        ok = sum(r["site_correct"] for r in rs)
        lig = sum(r["ligations"] for r in rs)
        print(f"  {arm:9s} ligations {lig}, template events {n}, complementary {ok}/{n} = {ok/max(n,1):.2f}")
        r2 = [r["round2"] for r in rs if "round2" in r]
        if r2:
            n2 = sum(x["site_events"] for x in r2); ok2 = sum(x["site_correct"] for x in r2)
            print(f"  {'round2':9s} templates {[x['template'] for x in r2]}, complementary {ok2}/{n2} = {ok2/max(n2,1):.2f}")


if __name__ == "__main__":
    main()
