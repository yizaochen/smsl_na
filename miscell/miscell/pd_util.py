import pandas as pd

def write_csv(f_csv, d_result, columns):
    df = pd.DataFrame(d_result)
    df = df[columns]
    df.to_csv(f_csv, index=False)
    print(f'Write DataFrame into {f_csv}')
