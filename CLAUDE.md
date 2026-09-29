# UK DRI scdownstream

## Use the repo agents for code changes

Implement pipeline changes through the two agents in `.claude/agents/`, not directly. This holds
for small changes too, and for follow-up fixes to code an agent wrote.

- **`pipeline-developer`**: `modules/local/**` (process bodies, Python/R templates,
  `environment.yml`), the Quarto report templates, and the `withName:` entries in
  `conf/modules.config`.
- **`nextflow-dsl2`**: `main.nf`, `workflows/*.nf`, `subworkflows/local/**`, the params in
  `nextflow.config`, `nextflow_schema.json` and `conf/test.config`.

When a change touches both, fix the interface first (process name, input order, emit names,
`task.ext.*` keys) and run the two agents in parallel. Afterwards, review and test their changes,
and update the docs (`README.md`, `docs/usage.md`, `docs/output.md`, `CHANGELOG.md`).

An agent that is resumed with a follow-up message keeps the instructions it started with. If
its definition file changed in the meantime, include the new rule in the message.
