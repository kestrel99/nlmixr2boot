# nlmixr2boot 0.2

* `runBootstrap()` now writes canonical shared `raw_results.*` artifacts, uses the common numbered run-directory and task-cache helpers from `nlmixr2utils`, and records a stable bootstrap master seed for restartable runs.

# nlmixr2boot 0.1

* Initial package split from `nlmixr2extra`, providing `runBootstrap()`,
  `sampling()`, `bootstrapNormalizedData()`, bootstrap summaries, tests, and the
  bootstrap vignette as a standalone package depending on `nlmixr2utils`.
