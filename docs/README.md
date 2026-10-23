# MethylSeq: Documentation

The documentation is published at [vibtobilab.github.io/MethylSeq](https://vibtobilab.github.io/MethylSeq/), built with [Astro](https://astro.build/) and [Starlight](https://starlight.astro.build/). Its pages are in [`src/content/docs/`](src/content/docs/):

- [Home](src/content/docs/index.mdx)
  - An overview of the pipeline and a quick start.
- [Usage](src/content/docs/usage.md)
  - How the pipeline works, how to run it and how to configure its optional steps.
- [Parameters](src/templates/parameters.md)
  - All the parameters of the pipeline, generated from `nextflow_schema.json` when the website is built.
- [Output](src/content/docs/output.md)
  - The results produced by the pipeline and how to interpret them.

The parameters, changelog and citations pages are generated from `nextflow_schema.json`, `CHANGELOG.md` and `CITATIONS.md` by [`src/integrations/generated-pages.mjs`](src/integrations/generated-pages.mjs). Images shown in the pages are in [`src/assets/`](src/assets/), files linked as they are (e.g. the interactive metro map) in [`public/`](public/).

To preview the website locally (Node.js >= 22), run from this directory:

```bash
npm ci
npm run dev
```

`npm run build` builds the website in `dist/` and fails on broken links and anchors.
