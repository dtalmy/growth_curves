import pandas as pd
import numpy as np
from sklearn import datasets, linear_model, metrics 


def linear_prediction(df):
    #df = pd.read_csv("drews_example_data.csv")
    xmin = df['x'].min()
    xmax = df['x'].max()
    x = np.array(df['x']).reshape((-1,1))
    y = np.array(df['y']).reshape((-1,1))
    reg = linear_model.LinearRegression() 
    reg.fit(x, y) 
    pred = reg.predict(np.linspace(xmin,xmax,len(df['x'])).reshape((-1,1)))
    df['prediction'] = pred.reshape(len(pred))
    return(df)
