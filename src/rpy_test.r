library("reticulate")

df = read.csv("drews_example_data.csv")

source_python("rpy_test.py")

df_python = linear_prediction(df)

print(df_python)
