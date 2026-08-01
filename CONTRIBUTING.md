# Contributing

Contributions are welcome. New inference backends should use `new_inference_model()` and return an `ecological_network`. Please add focused tests, document exported functions, run `devtools::check()`, and avoid coupling downstream metrics to a particular inference engine.

