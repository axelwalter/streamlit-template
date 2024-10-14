import streamlit as st
import re
from decimal import Decimal
import pandas as pd
import numpy as np

exact_masses = {'He': 4.002603,
                'Li': 7.016005,
                'Be': 9.012183,
                'B': 11.009305,
                'F': 18.998403,
                'Mg': 23.985042,
                'Al': 26.981541,
                'Si': 27.976928,
                'Cl': 34.968853,
                'Fe': 55.934942,
                'Cu': 62.929601,
                'Co': 58.933200,
                'Ni': 57.935348,
                'Zn': 63.929145,
                'Br': 78.918336,
                'C': 12.0,
                'O': 15.994915,
                'H': 1.007825,
                'N': 14.003074,
                'P': 30.973762,
                'S': 31.972071,
                'K': 38.963707,
                'Na': 22.989770,
                'Ca': 39.962591}


def get_element_dict(formula):
    elements = {}
    matches = re.findall(r'[A-Z][a-z]?[0-9]*', formula)
    if len("".join(matches)) != len(formula):
        return 'invalid formula'
    for m in matches:
        result = re.search(r'(\D+)(\d*)', m)
        if result[1] not in exact_masses:
            return 'invalid formula'
        if result[1] in elements:
            elements[result[1]] += int(result[2]) if result[2] != '' else 1
        else:
            elements[result[1]] = int(result[2]) if result[2] != '' else 1

    return elements


def get_formula_from_dict(elements):
    formula = ''
    for key, value in elements.items():
        if value == 1:
            formula += key
        else:
            formula += key + str(value)
    return formula


class Compound:
    def __init__(self, formula, name=None, charge=0, adducts={}):
        self.formula = formula
        self.elements = get_element_dict(formula)
        self.name = name
        self.charge = charge
        self.adducts = adducts

    def calc_mass(self, round_by=6):
        mass = Decimal('0')
        # calculate neutral mass
        for key, value in self.elements.items():
            mass += Decimal(str(exact_masses[key])) * value
        # add mass of adducts
        for key, value in self.adducts.items():
            mass += Decimal(str(exact_masses[key])) * value
        # remove electron masses if charge is positive
        if self.charge > 0:
            mass -= Decimal('0.000549') * self.charge
        # remove proton masses if charge is negative (in negative always remove protons)
        elif self.charge < 0:
            # remove H
            mass -= Decimal(str(exact_masses['H'])) * abs(self.charge)
            # and add back electrons
            mass += Decimal('0.000549') * abs(self.charge)
        if self.charge:
            # to get m/z, divide by charge
            mass /= abs(self.charge)
        return round(mass, round_by)

    def add_elements(self, add):
        to_add = get_element_dict(add)
        for key, value in to_add.items():
            if key in self.elements.keys():
                self.elements[key] += value
            else:
                self.elements[key] = value
            if self.elements[key] == 0:
                del self.elements[key]
        self.formula = get_formula_from_dict(self.elements)
        return self

    def del_elements(self, delete):
        to_delete = get_element_dict(delete)
        for key, value in to_delete.items():
            self.elements[key] -= value
            if self.elements[key] == 0:
                del self.elements[key]
        self.formula = get_formula_from_dict(self.elements)
        return self

    def check_formula(self):
        if self.elements == "invalid formula":
            return False
        if re.search(r'[0-9]+[a-z]+',self.formula) or re.search(r'^[a-z]+',self.formula) or re.search(r'[A-Z][a-z]{2,}',self.formula) or re.search(r'^[0-9]+',self.formula):
            return False
        for char in self.formula:
            if not char.isalnum():
                return False
        for key in self.elements.keys():
            try:
                exact_masses[key]
            except KeyError:
                return False
        return True
            
    def copy(self):
        return Compound(self.formula, self.name)
    
    def add_compound(self, compound, elimination = 'H2O', name='New Added Compound'):
        new_compound = self.copy()
        new_compound.name = name
        new_compound.add_elements(compound.formula)
        new_compound.del_elements(elimination)
        return new_compound
    
    def del_compound(self, compound, elimination = 'H2O', name='New Substracted Compound'):
        new_compound = self.copy()
        new_compound.name = name
        new_compound.del_elements(compound.formula)
        new_compound.add_elements(elimination)
        return new_compound
    
    def multiply(self, multiply, name='New Multiplied Compound', elimination = 'H2O'):
        new_compound = self.copy()
        new_compound.name = name
        new_formula = self.formula * multiply
        new_compound.elements = get_element_dict(new_formula)
        new_compound.formula = get_formula_from_dict(new_compound.elements)
        new_compound.del_elements(elimination*(multiply-1))
        return new_compound

    def change_name(self, name):
        self.name = name
        return self

def check_formula(formula):
    if re.search(r'[0-9]+[a-z]+',formula) or re.search(r'^[a-z]+',formula) or re.search(r'[A-Z][a-z]{2,}',formula) or re.search(r'^[0-9]+',formula):
        return False
    for char in formula:
        if not char.isalnum():
            return False
    for key in exact_masses.keys():
        try:
            check = exact_masses[key]
        except KeyError:
            return False
    return True


def get_mass(formula, adduct):
    if adduct == "[M+H]+":
        c = Compound(formula, charge=1)
        return c.calc_mass()
    if adduct == "[M+Na]+":
        c = Compound(formula, charge=1)
        try:
            c.add_elements("Na")
            c.del_elements("H")
            return c.calc_mass()
        except:
            return
    if adduct == "[M+2H]2+":
        c = Compound(formula, charge=2)
    if adduct == "[M-H2O+H]+":
        c = Compound(formula, charge=1)
        try:
            c.del_elements("H2O")
            return c.calc_mass()
        except:
            return
    if adduct == "[M-H]-":
        c = Compound(formula, charge=-1)
        return c.calc_mass()
    if adduct == "[M-2H]2-":
        c = Compound(formula, charge=-2)
        return c.calc_mass()
    if adduct == "[M-H2O-H]-":
        c = Compound(formula, charge=-1)
        try:
            c.del_elements("H2O")
            return c.calc_mass()
        except:
            return

def create_compound(formula, charge, adducts, neutral_loss, name, add_adduct_info):
    if not formula:
        st.session_state["mz_calc_error"] = "Enter sum formula for new metabolite."
        return pd.DataFrame()
    # try to build a valid compound with given input
    # check if adducts and charge are valid
    # get sum of adducts and compare with number of charges
    adducts = adducts.dropna()
    try:
        adducts["adduct"] = adducts["adduct"].astype(str)
        adducts["number"] = adducts["number"].astype(int)
    except:
        return pd.DataFrame()
    if adducts["number"].sum() > abs(charge) and charge != 0:
        st.session_state["mz_calc_error"] = f"Total number of adducts ({adducts['number'].sum()}) can not be greater then number of charges ({abs(charge)})."
        return pd.DataFrame()
    adduct_compound = Compound("", "adduct")
    if adducts["number"].sum() < charge and charge > 0:
        # fill adducts with protons up to charge number
        adduct_compound.add_elements(f"H{charge-adducts['number'].sum()}")
    # build adduct compound
    for adduct, number in zip(adducts["adduct"], adducts["number"]):
        c = Compound(adduct)
        if not c.check_formula():
            st.session_state["mz_calc_error"] = f"Invalid **adduct** sum formula: **{c.formula}**"
            return pd.DataFrame()
        for _ in range(int(number)):
            adduct_compound.add_elements(c.formula)
    # check if neutral losses are valid
    losses_compound = Compound(neutral_loss, "neutral loss")
    if not losses_compound.check_formula():
        st.warning(f"Invalid **neutral losses** sum formula: **{losses_compound.formula}**")
        return pd.DataFrame()
    # create compound here
    compound = Compound(formula, name, charge, adduct_compound.elements)
    # check formula
    if not compound.check_formula():
        st.session_state["mz_calc_error"] = f"Invalid sum formula: **{compound.formula}**"
        return pd.DataFrame()
    # check if there are enough protons to remove in negative mode
    if (("H" in compound.elements and compound.elements["H"] < abs(charge)) or ("H" not in compound.elements)) and charge < 0:
        st.session_state["mz_calc_error"] = f"Cannot remove protons from formula {compound.formula} in negative mode (negative charge)."
        return pd.DataFrame()
    # remove neutral losses
    compound.del_elements(losses_compound.formula)
    # adduct string
    # adduct string without charge
    if not charge:
        adduct_string = "[M]"
        if neutral_loss:
            adduct_string = adduct_string[:-1] + f"-{neutral_loss}]"
        if adduct_compound.formula:
            adduct_string = adduct_string[:-1] + f"+{adduct_compound.formula}]"
    else:
        # add all elements from adduct
        if adduct_compound.formula:
            adduct_notation = f"+{adduct_compound.formula}"
        else:
            adduct_notation = ""
        if charge == -1: # add substracted number of protons
            adduct_notation += f"-H"
        elif charge < -1:
            adduct_notation += f"-{abs(charge)}H"
        # determine charge sign
        if charge == 1:
            charge_sign = "+"
        elif charge > 1: 
            charge_sign = f"{charge}+"
        elif charge == -1:
            charge_sign = "-"
        elif charge < -1:
            charge_sign = f"{abs(charge)}-"
        if neutral_loss:
            adduct_string = f"[M-{neutral_loss}{adduct_notation}]{charge_sign}"
        else:
            adduct_string = f"[M{adduct_notation}]{charge_sign}"
    # determine a name
    # Determine compound name
    compound_name = name if name else f"{formula}#{adduct_string}"
    if name and add_adduct_info:
        compound_name += f"#{adduct_string}"
    # determine m/z
    mz = compound.calc_mass()

    # status.update(label="Valid metabolite created!", state="complete", expanded=True)
    return pd.DataFrame({"name": [compound_name], "sum formula": [formula], "adduct": [adduct_string], "mz": [mz], "RT": [np.nan], "peak width": [np.nan], "comment": [""]})

def can_eliminate(compound_one, compound_two, elimination_formula):
    # special case H2O, both compounds need to have OH groups
    if elimination_formula == "H2O":
        if not ("H" in compound_one and "O" in compound_one) or not ("H" in compound_two and "O" in compound_two):
            return False
    # Combine both compounds
    combined_compounds = {}
    # Combine the element counts from compound_one and compound_two
    for compound in [compound_one, compound_two]:
        for element, count in compound.items():
            combined_compounds[element] = combined_compounds.get(element, 0) + count

    # Check if there are at least as many of each element in elimination_compound
    for element, count in Compound(elimination_formula).elements.items():
        if element not in combined_compounds or combined_compounds[element] < count:
            return False
    
    return True

def build_compound(builder, charge, adducts, name, df, elimination, add_adduct_info, add_both_adducts):
    # check elimination
    if not Compound(elimination).check_formula():
        st.session_state["mz_calc_error"] = "Elimination product has invalid formula."
        return pd.DataFrame()
    # build compound
    total = {}
    for _, entry in builder.iterrows():
        formula = df[df["name"] == entry["metabolite"]]["sum formula"].iloc[0]
        if formula in total.keys():
            total[formula] += entry["number"]
        else:
            total[formula] = entry["number"]
    to_add = [list(c) for c in total.items() if c[1] > 0]
    to_sub = [list(c) for c in total.items() if c[1] < 0]
    # check if there is anything to add
    if not to_add:
        st.session_state["mz_calc_error"] = "Nothing to add."
        return pd.DataFrame()
    # start compound with first one to add
    compound = Compound(to_add[0][0])
    # reduce the number from first compound to add, already in starting compound
    to_add[0][1] -= 1
    # add the rest
    for entry in to_add:
        for _ in range(int(entry[1])):
            # check if it's possible to remove elemination product
            if not can_eliminate(compound.elements, Compound(entry[0]).elements, elimination):
                st.session_state["mz_calc_error"] = f"Can't add **{entry[0]}** to **{compound.formula}**. Impossible to remove elimination product **{elimination}**."
                return pd.DataFrame()
            compound = compound.add_compound(Compound(entry[0]), elimination=elimination)
    # remove all the substracted
    for entry in to_sub:
        for _ in range(abs(int(entry[1]))):
            compound = compound.del_compound(Compound(entry[0]), elimination=elimination)

    if not name:
        name = "+".join([f"{int(entry['number'])}({entry['metabolite']})" for _, entry in builder.iterrows()]).replace("+-", "-")

    if add_both_adducts and adducts["number"].sum() > 0:
        return pd.concat([create_compound(compound.formula, charge, pd.DataFrame({"adduct": [], "number": []}), "", name, True),
                        create_compound(compound.formula, charge, adducts, "", name, True)])
    return create_compound(compound.formula, charge, adducts, "", name, add_adduct_info)


def save_df(new_compound_df, input_table_path): 
    df = pd.read_csv(input_table_path)
    if not new_compound_df.empty:
        if df["name"].isin([new_compound_df["name"][0]]).any():
            st.session_state["mz_calc_error"] = "Metabolite names have to be unique."
            if "mz_calc_success" in st.session_state:
                del st.session_state["mz_calc_success"]
            st.rerun()
        pd.concat([new_compound_df, df]).to_csv(input_table_path, index=False)
        st.session_state["mz_calc_success"] = [f"**{row['name']}** with adducts **{row['adduct']}** and m/z **{row['mz']}**" for _, row in new_compound_df.iterrows()]
        if "mz_calc_error" in st.session_state:
            del st.session_state["mz_calc_error"]
        st.rerun()
    else:
        if "mz_calc_success" in st.session_state:
            del st.session_state["mz_calc_success"]
        st.rerun()

def validate_dataframe(df):
    # Define expected dtypes and column names
    expected_dtypes = {
        "name": str, 
        "sum formula": str, 
        "adduct": str, 
        "mz": float, 
        "RT": float, 
        "peak width": float, 
        "comment": str
    }
    # Try to convert data types
    for col, dtype in expected_dtypes.items():
        try:
            df[col] = df[col].astype(dtype)
        except ValueError:
            st.session_state["mz_calc_error"] = f"Error: Column '{col}' cannot be converted to {dtype.__name__}"
            return False
        except KeyError:
            st.session_state["mz_calc_error"] = f"Error: Column '{col}' is missing in the uploaded file."
            return False
    return True


HELP = """
The m/z calculator facilitates the calculation of mass-to-charge ratios (m/z) for metabolites

#### ➕ New Metabolite Tab
In this section, you can add new metabolites to your analysis.

- **Sum Formula**: Enter the molecular formula of the new metabolite.
- **Metabolite Name**: Optionally provide a name for the metabolite. If left blank, a name will be generated automatically.
- **Neutral Losses**: Optionally specify any neutral losses in the sum formula format (e.g., H2O).
- **Charge**: Input the charge of the metabolite. Use negative numbers for negative ion mode and positive numbers for positive ion mode.
- **Adducts**: Define adducts (excluding protons) that correspond to the number of charges. In positive mode, any remaining charges will be considered as protons. In negative mode, a number of protons equal to the absolute charge number will be removed, irrespective of additional specified adducts.

After filling out the form, click "Add new metabolite" to calculate the m/z and add the metabolite to your table.

#### 📟 Combine Metabolites Tab
Use this section to combine existing metabolites into larger molecules.

- **Metabolites to Combine**: Select from pre-existing metabolites in your table to combine them. Specify the number of times each metabolite is to be used in the combination.
- **Charge**: Set the overall charge for the combined metabolite.
- **Metabolite Name**: Optionally provide a name for the combined metabolite. If left blank, a name will be generated automatically.
- **Elimination Product**: Optionally specify any elimination product to be removed when combining metabolites (default is "H2O").
- **Adducts**: Like in the "New Metabolite" tab, define adducts for the combined metabolite.

#### ⬆️ Upload existing table Tab

- upload a m/z table which has been generated with this tool (e.g. to get a table from another workspace)

After setting up your combination, click "Calculate metabolite" to process and add it to your table.

## Editing and Saving Your Work
Below the tabs, you'll find an editable table displaying your metabolites. You can modify this table as needed. Remember that metabolite names need to be unique. After making changes, the tool will automatically save your updated data.

"""