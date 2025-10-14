import pandas as pd
import logging
pd.options.display.max_columns = None
pd.options.display.max_rows = None

logger = logging.getLogger(__name__)

def read_excel(infile,sheet_name=None):
    '''
    For reading excel as pd.DataFrame
    '''
    if sheet_name is None:
        df_models = pd.read_excel(infile)
    else:
        df_models = pd.read_excel(infile, sheet_name= sheet_name, engine='openpyxl')
    return df_models

def write_excel(infile,df_models):
    '''
    For writing excel as pd.DataFrame
    '''    
    df_models.to_excel(infile, sheet_name='Sheet1', index=False)

def excel_sheet_extraction(infile,sheet_name):
    '''
    For extraction of one sheet from execl file.
    '''
    excel_file = pd.ExcelFile(infile)
    df_models = excel_file.parse(sheet_name)  
    return df_models

def excel_sheet_replace(infile,sheet_name,df_models):
    '''
    For replacing the sheet in excel
    '''
    with pd.ExcelWriter(infile, engine='openpyxl', \
                        mode='a', if_sheet_exists='replace') as writer:
        df_models.to_excel(writer, sheet_name=sheet_name, index=False)




