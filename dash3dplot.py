from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import pandas as pd

app = Dash(__name__)

cycleList = [i for i in range(10000,20000,1000)]

app.layout = html.Div([
    html.Div(
        children=[
            html.H3("Cycle Number",style={"text-align":"center"}),
            html.Div([
                dcc.Dropdown(cycleList, min(cycleList), id='cycle-dropdown',clearable=False),
            ],style={"width":"50%","margin": "0 auto"}),
            dcc.Graph(id='graph-with-dropdown',style={"text-align": "center","height":"90vh"})
        ],style={'display': 'inline-block',"width":"48vw"} 
    ),
    html.Div(
        children=[
            html.H3("Projection direction",style={"text-align":"center"}),
            html.Div([
                dcc.Dropdown(['X','Y','Z'], 'Z', id='cycle-dropdown-2',clearable=False)
            ],style={"width":"50%","margin": "0 auto"}),
            dcc.Graph(id='graph-with-dropdown-2',style={"text-align": "center","height":"90vh"})
        ],style={'display': 'inline-block',"width":"48vw"} 
    )
])

@app.callback(
    Output("graph-with-dropdown", "figure"), 
    Input("cycle-dropdown", "value"))
def update_bar_chart(selected_cycle):
    # df = pd.read_csv(f'./positions/position_-5000.00_{selected_cycle}.csv') # replace with your own data source
    df = pd.read_csv(f'./GrapheneResults/position_{selected_cycle}_1000000.00.csv') # replace with your own data source
    df.columns = ['x','y','z']
    fig = px.scatter_3d(df, 
        x='x', y='y', z='z')
    # fig = px.scatter(df,x='x',y='y')
    return fig

@app.callback(
    Output("graph-with-dropdown-2", "figure"),
    Input("cycle-dropdown", "value"), 
    Input("cycle-dropdown-2", "value"))
def update_bar_chart_2(selected_cycle,projection_dir):
    # df = pd.read_csv(f'./positions/position_-5000.00_{selected_cycle}.csv') # replace with your own data source
    df = pd.read_csv(f'./GrapheneResults/position_{selected_cycle}_1000000.00.csv') # replace with your own data source
    df.columns = ['x','y','z']
    # fig = px.scatter_3d(df, 
    #     x='x', y='y', z='z')
    if projection_dir == 'Z':
        fig = px.scatter(df, x='x',y='y')
    elif projection_dir == 'Y':
        fig = px.scatter(df, x='x',y='z')
    else:
        fig = px.scatter(df, x='y',y='z')
    return fig

app.run_server(debug=False)