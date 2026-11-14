## What and why

<!-- One or two sentences: what this PR changes and which plate or part of the engine it touches. -->

## Checklist

Details for each line are in [CONTRIBUTING.md](../CONTRIBUTING.md#before-you-open-a-pull-request).

- [ ] Jury-facing pages start in English; every map label has an `en` version; internal Ukrainian text is in `docs/uk/`
- [ ] The plate claims no more than the tests show; "partly" and "no" results are said on the plate
- [ ] Each number has the same value everywhere it appears
- [ ] Tests behind claims are pre-registered in `notes/experiments.md`; after-the-fact choices are marked post-hoc
- [ ] New datasets are in `docs/data-guide.md` and the `README.md` credits, with licences
- [ ] A teammate can rebuild this from the repo (`sars fetch` or inputs described in the data guide)
- [ ] Commits use a GitHub noreply e-mail
- [ ] `uv run python tests/smoke.py` passes
