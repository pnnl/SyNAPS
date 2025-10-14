import pandas as pd

def update_config(config_df, updates):
    """
    Update specific settings in the configuration DataFrame.

    Parameters:
        config_df (pd.DataFrame): The default configuration DataFrame.
        updates (dict): Key-value pairs where key is column name and value is desired status ('sym'/'num'/'user').
    
    Returns:
        pd.DataFrame: Updated configuration DataFrame.
    """
    updated_df = config_df.copy()  # Make a copy to avoid modifying original DataFrame
    
    for key, value in updates.items():
        if key in updated_df.columns:
            updated_df.at[0, key] = value  # Update specific value for the first row
        else:
            print(f"Warning: {key} not found in configuration!")

    return updated_df