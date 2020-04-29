# importing copy module
import plotly.graph_objects as go
from plotly.subplots import make_subplots


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

    @staticmethod
    def GenerateBinaryCommunityHeatMap(Title, X_AxisTitle, Y_AxisTitle, TwoDimMatrix, PerimeterValues):
        # Cluster Heatmap plotting individual:
        pValues = PerimeterValues
        matrix = TwoDimMatrix
        fig = make_subplots(rows=1, cols=1, vertical_spacing=0.2, shared_xaxes=False, specs=[[{}]])
        fig.add_trace(
            go.Heatmap(y=pValues, x=pValues, z=matrix, colorscale="Spectral", reversescale=True, showscale=False), 1, 1)
        fig.update_layout(yaxis_title='Windows', title_text="<b>" + Title + "</b>")
        fig.show()