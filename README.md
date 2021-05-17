# mg-evaluation

The code in this repository demonstrates the Multi-Generation Evaluation approach documented in [Holistic Evaluation of Biodegradation Pathway Prediction: Assessing Multi-Step Reactions and Intermediate Products](https://chemrxiv.org/articles/preprint/Holistic_Evaluation_of_Biodegradation_Pathway_Prediction_Assessing_Multi-Step_Reactions_and_Intermediate_Products/14315963).

The repository also contains two files, `TRAIN_SOIL_PW.csv` and `TEST_SOIL_PW.csv`. They each contain a list of pathway names included respectively in the `TRAIN_SOIL` and `TEST_SOIL` packages, which were used in the experiment. 

## Installation

The code is intended to be used with Python 3.5.3 or later versions. Other required libraries are listed in `requirements.txt` and can be installed using the package manager [pip](https://pip.pypa.io/en/stable/).

```bash
pip install -r requirements.txt
```

## Usage

Simply run this in the terminal

```bash
python main.py
```
after parameters such as `obsPathwayQuery` and `predPathwayQuery` are configured at the bottom of the code script