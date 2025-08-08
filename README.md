# lab-data-processing

This project provides software to process and analyze data from diagnostic ELISA tests.

---

## Project Structure
```
lab-data-processing/
├── reference_match.py # Allows you to find the instrument/peer group that most closely matches your sample data.
├── sample_data.txt # Your ELISA sample values (input)
├── reference_data.txt # Reference instrument values (input)
├── requirements.txt # Python dependencies
└── README.md # Project description and usage
```

## Input Format

### `reference_data.txt`

A whitespace-separated table with one row per instrument per sample. ex:
```
Instrument Sample Number Group # Labs Mean SD Low Range High Range Uncertainty
Abbott Alinity i/Alinity reagent 1 IA-01 43 315.3 17.9 220 410 3.41
```

### `sample_data.txt`

A plain text file containing a single sample value per line, in order of sample number. ex:
```
115
223
137.1
893
20.5
```

## Installation & Running

1. **Clone this repository**:
```bash
git clone https://github.com/nataliekozlowskii/lab-data-processing.git
cd lab-data-processing
```
2. **Install dependencies:**
```
pip install -r requirements.txt
```
3. **Run the script:**
```
python reference_match.py
```
