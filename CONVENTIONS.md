# Repository conventions — gabrielgrant/monstertruck

This repository is a single GitHub fork serving **two** upstreams. GitHub permits one
fork per network, and `virtualritz/monstertruck` is itself a GitHub fork of
`ricosjp/truck` (this repo's API metadata: parent = virtualritz/monstertruck,
source = ricosjp/truck) — so this one fork can open pull requests against **both**
upstreams. This branch (`meta/conventions`) is documentation only: pushed to the fork,
never merged anywhere.

## Remotes (local clones)

| remote | URL | push? |
|---|---|---|
| `origin` | github.com/gabrielgrant/monstertruck | yes — the only push target |
| `monstertruck` | github.com/virtualritz/monstertruck | never push |
| `truck` | github.com/ricosjp/truck | never push |

## Branch namespaces

| namespace | based on | purpose |
|---|---|---|
| `mirror/monstertruck-master`, `mirror/truck-master` | the respective upstream master | pushed snapshots of the upstream masters; refresh with the commands below |
| `mt/<fix|feat>/<name>` | `monstertruck/master` | one self-contained change, PR-ready for **virtualritz/monstertruck** |
| `tr/<fix|feat>/<name>` | `truck/master` | one self-contained change, PR-ready for **ricosjp/truck** |
| `patches/monstertruck` | `monstertruck/master` | integration branch: master + all `mt/*` merged (each `--no-ff`) + the `.blueprints` submodule-removal commit (see caveat). Consumed by downstream (openshape) as a cargo git dependency. |
| `patches/truck` | `truck/master` | integration branch: master + all `tr/*` merged |
| `meta/*` | orphan | repo documentation, never merged |

Namespace prefixes are `mt/`/`tr/` (not `monstertruck/`/`truck/`) to avoid ref-shorthand
ambiguity with the remote names.

## Refresh procedure

```sh
git fetch monstertruck && git fetch truck
git push origin refs/remotes/monstertruck/master:refs/heads/mirror/monstertruck-master
git push origin refs/remotes/truck/master:refs/heads/mirror/truck-master
# rebuild an integration branch after upstream moves (example: monstertruck):
git checkout -B patches/monstertruck monstertruck/master
git rm --cached .blueprints && git config -f .gitmodules --remove-section submodule..blueprints \
  && git commit -m "chore: drop private .blueprints submodule reference (unfetchable by cargo git deps)"
for b in $(git branch --list 'mt/*' --format='%(refname:short)'); do git merge --no-ff "$b"; done
git push -f origin patches/monstertruck
```

## Caveats

1. **`.blueprints` submodule**: upstream monstertruck's `.gitmodules` references the
   PRIVATE `virtualritz/blueprints` repo. Cargo git-dependencies fetch every submodule
   unconditionally, so any branch consumed by cargo MUST carry the removal commit
   (`patches/monstertruck` does). Re-apply it on every rebase of that branch.
2. **PR targeting**: when opening a PR from an `mt/*` branch, set the base repository to
   virtualritz/monstertruck; from a `tr/*` branch, ricosjp/truck. GitHub defaults the
   base to the parent (virtualritz) — change it manually for `tr/*` PRs.
3. Cross-upstream duplicates: features developed on `mt/*` that also apply to truck
   must be hand-ported to a `tr/*` branch (the projects share ancestry but have
   diverged in naming/APIs); note the sibling branch in both PR bodies.
4. Considering renaming this GitHub repo to reflect its dual role: GitHub redirects
   old URLs (including git fetch/clone) after a rename, so downstream cargo git deps
   keep working, but update them promptly anyway.

## Current inventory (2026-07-09)

- `mt/fix/fillet-adjacent-edges`, `mt/fix/wasm-safe-clock`,
  `mt/fix/unit-circle-arc-accuracy`, `mt/feat/path-sweep`, `mt/feat/shell-thicken`,
  `mt/feat/draft`, `mt/feat/planar-shared-and-curved` (this last one bases on
  `patches/monstertruck`, not master — it refactors shell+draft; upstream it only
  after those land)
- `tr/fix/coincident-plane-booleans`
- `patches/monstertruck` (consumed by openshape), `patches/truck`
