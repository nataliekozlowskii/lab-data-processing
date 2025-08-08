import os
import pandas as pd
from scipy.spatial.distance import euclidean
import numpy as np

class DataComparison:
    """
    Compares experimental sample data with reference dataset of instrument measurements.

    This class offers methods to find the closest instrument match by:
        - Euclidean distance
        - Max number of samples falling within reference bounds
        - Max number of samples within a percent deviation from the reference mean

    Args:
        reference_df (pd.DataFrame): Reference dataframe with rows as instrument data 
                                     for a specific sample.
        sample_values_dict (dict): Dictionary mapping sample numbers to concentration values.
    """

    def __init__(self, reference_df, sample_values_dict):

        self.reference_df = reference_df
        self.sample_values_dict = sample_values_dict

    def find_closest_euclidean_distance(self):

        """
        Finds instrument from the reference dataset whose mean sample values
        have the smallest Euclidean distance from the provided sample data.

        Returns:
            tuple: (best_match, best_distance) where
                - best_match (str): Instrument name with the closest Euclidean distance.
                - best_distance (float): Minimum Euclidean distance.
        """

        sample_vector = np.array([self.sample_values_dict[i] for i in sorted(self.sample_values_dict)])

        # create a table with each row as an instrument,
        # each col as a sample number,
        # and each value as the mean for that instrument and sample
        pivot = self.reference_df.pivot_table(
            index="Instrument", 
            columns="Sample Number", 
            values="Mean"
        )

        # drop rows that don't have all samples from consideration
        pivot = pivot.dropna()

        # dict of {instrument: distance from sample vector}
        distances = {
            instrument: euclidean(pivot.loc[instrument].values, sample_vector)
            for instrument in pivot.index
        }

        best_match = min(distances, key=distances.get)
        return (best_match, distances[best_match])
    
    def count_within_bounds(sample_values, low_bounds, high_bounds):

        """Return how many of the sample values are within the given bounds."""

        return sum(low <= y <= high for y, low, high in zip(sample_values, low_bounds, high_bounds))
    
    def find_most_within_bounds(self):

        """
            Finds instrument from the reference dataset for which the greatest
            number of sample concentrations are within the reference range.

            Returns:
                tuple: (best_match, max_within_bounds) where
                    - best_match (str): Instrument name with the most within range samples.
                    - max_within_bounds (int): Highest number of within range samples.
        """

        instruments = self.reference_df["Instrument"].unique()
        best_match = None
        max_within_bounds = -1

        for instrument in instruments:
            # filter rows for this instrument
            df_instrument = self.reference_df[self.reference_df["Instrument"] == instrument]

            count = 0
            for _, row in df_instrument.iterrows():
                sample_num = row["Sample Number"]
                value = self.sample_values_dict.get(sample_num, None)

                if value is not None:
                    if row["Low Range"] <= value <= row["High Range"]:
                        count += 1

            if count > max_within_bounds:
                max_within_bounds = count
                best_match = instrument
        
        return best_match, max_within_bounds
    
    def find_most_within_percent(self, percent):

        """
            Finds instrument from the reference dataset for which the greatest
            number of sample concentrations are within a percent value of the mean
            reference value.

            Args:
                percent (float): Acceptable percent deviation from the mean, in percents.

            Returns:
                tuple: (best_match, max_within_bounds) where
                    - best_match (str): Instrument name with the most within range samples.
                    - max_within_percent (int): Highest number of samples within a percent of the mean.
        """

        instruments = self.reference_df["Instrument"].unique()
        best_match = None
        max_within_percent = -1

        for instrument in instruments:
            # filter rows for this instrument
            df_instrument = self.reference_df[self.reference_df["Instrument"] == instrument]

            count = 0
            for _, row in df_instrument.iterrows():
                sample_num = row["Sample Number"]
                value = self.sample_values_dict.get(sample_num, None)

                if value is not None:
                    percent_diff = (value - row["Mean"])/row["Mean"]
                    if abs(percent_diff) <= percent/100:
                        count += 1

            if count > max_within_percent:
                max_within_percent = count
                best_match = instrument
        
        return best_match, max_within_percent
    
class ReferenceData:
    """
    Loads reference data into a pd.DataFrame and provides methods to access reference data.

    Args:
        reference_data_file_path (string): Path to the reference data file.
    """

    def __init__(self, reference_data_file_path):
        # {Sample Name: reading that section of the data table (True/False)}
        self.in_sample_dict = {}
        for sample in range(1, 11):
            self.in_sample_dict[sample] = False

        # {Group Name: reading that section of the data table (True/False)}
        self.in_group_dict = {}
        for group in ["Peer Group", "Instrument Group", "Method Group",  "Reagent Group"]:
            self.in_group_dict[group] = False

        self.reference_data_file_path = reference_data_file_path
        self.reference_df = pd.DataFrame(self.load_reference_data())
    
    def set_in_sample(self, sample):

        """
        Specify the sample number whose data we are currently 
        reading in the reference table.
        """

        self.in_sample_dict[sample] = True
        for i in range(1, 11):
            if i != sample:
                self.in_sample_dict[i] = False
    
    def get_in_sample(self):

        """
        Return the sample number whose data we are currently 
        reading in the reference table.
        """

        for sample in self.in_sample_dict.keys():
            if self.in_sample_dict[sample] == True:
                return sample

        return 0

    def set_in_group(self, group):

        """
        Specify the group type section whose data we are currently 
        reading in the reference table.
        """

        self.in_group_dict[group] = True
        for group_key in self.in_group_dict.keys():
            if group_key != group:
                self.in_group_dict[group_key] = False

    def get_in_group(self):

        """
        Return the group type section whose data we are currently 
        reading in the reference table.
        """

        for group in self.in_group_dict.keys():
            if self.in_group_dict[group] == True:
                return group
            
        return "No Group Found"

    def get_reference_df(self):
        return self.reference_df.copy()

    def load_reference_data(self):
        """
        Creates a pd.DataFrame of the reference data 
        where each row is a new instrument/sample match.

        Args:
            None

        Returns:
            reference_df (pd.DataFrame): Dataframe containing reference data with columns:
                - "Instrument" (str) --> name of the instrument
                - "Sample Number" (int64) --> sample number 1-10
                - "Group" (str) --> "Peer Group", "Instrument Group"
                            "Method Group", or "Reagent Group" (str)
                - "# Labs" (int64) --> # of labs used in data collection (int64)
                - "Mean" (float64) --> mean of the data
                - "SD" (float64) --> standard deviation of the data
                - "Low Range" (int64) --> lower bound for range of the data
                - "High Range" (int64) --> upper bound for range of the data
                - "Uncertainty" (float64) --> standard uncertainty of the measurements
        """

        with open(self.reference_data_file_path, "r") as f:
            all_lines = f.readlines()

        reference_data_list = []

        for line in all_lines:
            if "SAMPLE IA-01" in line:
                self.set_in_sample(1)
            elif "SAMPLE IA-02" in line:
                self.set_in_sample(2)
            elif "SAMPLE IA-03" in line:
                self.set_in_sample(3)
            elif "SAMPLE IA-04" in line:
                self.set_in_sample(4)
            elif "SAMPLE IA-05" in line:
                self.set_in_sample(5)
            elif "SAMPLE IA-06" in line:
                self.set_in_sample(6)
            elif "SAMPLE IA-07" in line:
                self.set_in_sample(7)
            elif "SAMPLE IA-08" in line:
                self.set_in_sample(8)
            elif "SAMPLE IA-09" in line:
                self.set_in_sample(9)
            elif "SAMPLE IA-10" in line:
                self.set_in_sample(10)
            elif "Peer Group" in line:
                self.set_in_group("Peer Group")
            elif "Instrument Groups" in line:
                self.set_in_group("Instrument Group")
            elif "Method Groups" in line:
                self.set_in_group("Method Group")
            elif "Reagent Groups" in line:
                self.set_in_group("Method Group")
            elif line.strip() == "" or "All Participants" in line:
                continue
            else: # we are in a data block
                sample_num = self.get_in_sample()
                if sample_num == 0:
                    print("Error. Could not identify the sample block number")
                else:
                    # remove the hyphen indicating range
                    line = line.replace("-", "")
                    line_data = line.split()
                    instrument_name = " ".join(line_data[:-6])
                    sample_num = self.get_in_sample()
                    group = self.get_in_group()
                    row_data = [instrument_name, sample_num, group] + line_data[-6:]
                    reference_data_list.append(row_data)

        reference_df = pd.DataFrame(reference_data_list)
        reference_df.columns = ["Instrument", 
                     "Sample Number", 
                     "Group", 
                     "# Labs", 
                     "Mean", 
                     "SD", 
                     "Low Range", 
                     "High Range", 
                     "Uncertainty"]
        
        # convert numeric columns to numeric types
        numeric_cols = ["Sample Number", "# Labs", "Mean", "SD", "Low Range", "High Range", "Uncertainty"]
        reference_df[numeric_cols] = reference_df[numeric_cols].apply(pd.to_numeric, errors='coerce')

        return reference_df

    def get_instrument_values(self, instrument):
        return self.reference_df[self.reference_df["Instrument"] == instrument]
    
class SampleData:
    """
    Loads sample data into a dictionary and provides methods to access sample data.
    
    Args:
        sample_data_file_path (string): Path to the sample data file.
    """

    def __init__(self, sample_data_file_path):

        self.sample_data_file_path = sample_data_file_path
        self.sample_values_dict = self.load_sample_data()
    
    def load_sample_data(self):

        """Load sample concentrations into a dictionary keyed by sample number."""

        sample_values_dict = {}
        with open(self.sample_data_file_path, "r") as f:
            all_lines = f.readlines()
        
        i = 1
        for line in all_lines:
            sample_values_dict[i] = float(line.strip())
            i += 1
        
        return sample_values_dict

    def get_sample_values_dict(self):
        return self.sample_values_dict.copy()
            
def main(reference_data_file_name, sample_data_file_name, within_percent):

    reference_data_file_path = os.path.join(os.getcwd(), reference_data_file_name)
    sample_data_file_path = os.path.join(os.getcwd(), sample_data_file_name)

    reference_loader = ReferenceData(reference_data_file_path)
    reference_df = reference_loader.get_reference_df()

    sample_loader = SampleData(sample_data_file_path)
    sample_values_dict = sample_loader.get_sample_values_dict()

    data_compare = DataComparison(reference_df, sample_values_dict)

    best_euclidean_instrument, best_euclidean_distance = data_compare.find_closest_euclidean_distance()
    print(f"Instrument with lowest Euclidean distance: {best_euclidean_instrument} "
          f" with distance {best_euclidean_distance:.4f}")

    best_within_bounds_instrument, count_within_bounds = data_compare.find_most_within_bounds()
    print(f"Instrument with most sample values within reference range: "
          f"{best_within_bounds_instrument} with max within {count_within_bounds:.4f}")

    best_within_percent_instrument, count_within_percent = data_compare.find_most_within_percent(within_percent)
    print(f"Instrument with most sample values within {within_percent}: "
          f" {best_within_percent_instrument} with most {count_within_percent:.4f}")


if __name__ == "__main__":

    reference_data_file_name = "reference_data.txt"
    sample_data_file_name = "sample_data.txt"
    within_percent = 30

    main(reference_data_file_name, sample_data_file_name, within_percent)