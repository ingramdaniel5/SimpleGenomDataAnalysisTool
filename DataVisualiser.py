from Profile import GenomeProfile as Profile
from Window import GenomeWindow as Window
from WindowDataSet import WindowDataSet as WSet
from HeatMapWindow import heatmap
import matplotlib.pyplot as plt
import random
# importing copy module
import copy
import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px


class DataVisualiser:

    @staticmethod
    def GenerateHeatMap(Title, X_AxisTitle, Y_AxisTitle, TwoDimMatrix, PerimeterValues):
        # Cluster Heatmap plotting individual:
        pValues = PerimeterValues
        matrix = TwoDimMatrix
        fig = make_subplots(rows=1, cols=1, vertical_spacing=0.2, shared_xaxes=False, specs=[[{}]])
        fig.add_trace(go.Heatmap(y=pValues, x=pValues, z=matrix, colorscale="Spectral", reversescale=True, showscale=False), 1, 1)
        fig.update_layout(yaxis_title='Radial score', title_text="<b>" + Title + "</b>")
        fig.show()