# [1.2.1] - 2021-06-29
* **Fixed** an issue (#14) where bcSignatures would fail when used with pathways
# [1.2.0] - 2021-06-8
* **Fixed** an issue where bcMerge would fail if the background matrix was empty. 
# [1.1.1] - 2021-05-20
* **Updated** pathways with the latest sources from MSigDB 7.3.
* **Added** a directory called `data-raw` with both the code and the data used to generate the pathway signatures. 

# [1.1] - 2021-04-14

* **Added** two new arguments in `bcUMAP`: `method` and `return.model`. `method` can be `"umap-learn"` or `"uwot"` (default). `return.model` is a logical argument that determines whether the function will return the uwot model (only valid when `method = "uwot"`).
* **Changed tutorials**: Typos have been corrected, chunks of code now specify the UMAP method (`"umap-learn"`) and figures have been redrawn
(they look the same as the originals).

# [1.0] - 2021-04-8

* Initial public release: subsequent Beyondcell releases will follow [Semantic Versioning guidelines](https://semver.org/).
