# MethylSeq: Documentation

The documentation is published at [vibtobilab.github.io/MethylSeq](https://vibtobilab.github.io/MethylSeq/) and is split into the following pages:

- [Home](index.md)
  - An overview of the pipeline and a quick start.
- [Usage](usage.md)
  - How the pipeline works, how to run it and the modules added to nf-core/methylseq.
- [Parameters](parameters.md)
  - All the parameters of the pipeline, generated from `nextflow_schema.json` when the website is built.
- [Output](output.md)
  - The results produced by the pipeline and how to interpret them.

To preview the website locally, run from the root of the repository:

```bash
pip install -r docs/requirements.txt
mkdocs serve
```
