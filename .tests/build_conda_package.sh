
if [ -z "$1" ]
then
    PYTHON_VERSION=3.6
else
    PYTHON_VERSION=$1
fi

git clone https://gitlab.com/bu_cnio/beyondcell_conda_recipe
conda mambabuild beyondcell_conda_recipe/r-beyondcell --output-folder ./
mamba install --use-local --update-deps r-beyondcell
