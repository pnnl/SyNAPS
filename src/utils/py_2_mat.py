# from sympy import symbols, Matrix
import textwrap

# Convert expression to MATLAB syntax
def sympy_to_matlab(expr):
    return str(expr).replace("**", "^").replace("*", ".*")

# Format matrix to multiline MATLAB string with optional line wrapping
def matrix_to_matlab(matrix, var_name='A', line_limit=80):
    matlab_rows = []
    for row in matrix.tolist():
        row_exprs = [sympy_to_matlab(e) for e in row]
        row_str = ', '.join(row_exprs)

        # Wrap long lines with continuation (...) in MATLAB
        if len(row_str) > line_limit:
            wrapped = textwrap.wrap(row_str, width=line_limit, break_long_words=False)
            wrapped = ["    " + w + " ..." for w in wrapped[:-1]] + ["    " + wrapped[-1]]
            matlab_rows.append("\n".join(wrapped))
        else:
            matlab_rows.append("    " + row_str)
    
    body = ";\n".join(matlab_rows)
    return f"{var_name} = [\n{body}\n];"

# Optional: add symbolic declarations
def matlab_syms(vars):
    return "syms " + ' '.join(str(v) for v in vars) + "\n"

def matlab_value_assignments(var_values, var_name='A'):
    # Convert symbolic vars to strings
    assigns = [f"{str(var)}_val = {val};" for var, val in var_values.items()]
    
    var_list = [str(var) for var in var_values.keys()]
    val_list = [f"{str(var)}_val" for var in var_values.keys()]
    
    subs_line = f"{var_name}_sub = subs({var_name}, [{', '.join(var_list)}], [{', '.join(val_list)}]);"
    numeric_line = f"{var_name}_numeric = double({var_name}_sub);"
    
    return "\n".join(assigns + ["", subs_line, numeric_line])