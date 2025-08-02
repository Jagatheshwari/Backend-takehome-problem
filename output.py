[tool.poetry]
name = "get-papers"
version = "0.1.0"
description = "Fetch PubMed papers with pharma affiliations"
authors = ["Your Name <you@example.com>"]

[tool.poetry.dependencies]
python = "^3.8"
biopython = "^1.81"

[tool.poetry.scripts]
get-papers-list = "cli:main"

[build-system]
requires = ["poetry-core"]
build-backend = "poetry.core.masonry.api"
