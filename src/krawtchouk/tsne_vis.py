import os
import json
import numpy as np
import pandas as pd
import cv2
import glob
from sklearn.manifold import TSNE
from krawtchouk_moments import krawtchouk_moments

import base64
import dash
from dash import Dash, dcc, html, Input, Output, no_update
import plotly.express as px
import plotly.graph_objects as go

def read_data(data_dir, N, M):

    images = np.zeros((len(glob.glob(data_dir + '/*/*/*.png')), N, M))
    labels = np.zeros(len(images), dtype='int')
    chars = [None] * len(glob.glob(data_dir + '/*/'))

    i = 0
    for label, folder in enumerate(glob.glob(data_dir + '/*/')):
        tmp = folder.split('/')
        char = tmp[-2]
        chars[label] = char
        for file in glob.glob(data_dir + '/{}/*/*.png'.format(char)):
            img = cv2.resize(cv2.imread(file), (M, N))[:, :, 0]
            images[i, :, :] = img
            labels[i] = label
            i += 1

    return images, labels, chars


def extract_krawtchouk_moments(images, N_order, M_order):

    p = [0.5, 0.5]
    N, M = images[0].shape
    km = np.zeros(images.shape)
    for i in range(len(images)):
        km[i, :, :] = krawtchouk_moments(images[i], N, M, p)

    km = km[:, :N_order, :M_order]
    km = km.reshape((len(km), N_order * M_order))
    return km


def get_image_url(img):
        _, buffer_img = cv2.imencode('.jpg', img)
        img_base64 = base64.b64encode(buffer_img).decode("utf-8")
        return "data:image/jpeg;base64, " + img_base64


def nns_to_elems(nns, dists, ranks):
    elems = []
    text_style = {'margin': '0px'}
    for i, nn in enumerate(nns):
        elems.append(html.Div([
            html.Img(
                src=get_image_url(images_masked[nn]),
                style={"width": "70px"},
            ),
            html.Div([
                html.P(
                    '{}. {} ({})'.format(str(ranks[i]), chars[int(df['label'][nn])], str(nn)),
                    style=text_style
                ),
                html.P(
                    'Dist: {}'.format(str(round(dists[nn], 2))),
                    style=text_style
                )
            ])
            
        ], style={
            'display': 'inline-block', 
            'white-space': 'normal',
            'width': '120px'
        }))
    return elems


if os.path.isdir('vis_cache'):
    with open('vis_cache/images.npy', 'rb') as f:
        images = np.load(f)
    with open('vis_cache/km.npy', 'rb') as f:
        km = np.load(f)
    with open('vis_cache/meta.json', 'r') as f:
        meta = json.load(f)
        N, M, N_order, M_order = meta['N'], meta['M'], meta['N_order'], meta['M_order']
        labels, chars = np.array([int(i) for i in meta['labels']]), meta['chars']
else:
    data_dir = '/Users/kx/Docs/github/git_cccr_sl/chinese-calligraphy-recognition-sl/data/cccr/train'
    N, M = 128, 128
    N_order, M_order = 20, 20
    images, labels, chars = read_data(data_dir, N, M)
    km = extract_krawtchouk_moments(images, N_order, M_order)
    meta = {
        'N': N,
        'M': M,
        'N_order': N_order,
        'M_order': M_order,
        'labels': [str(i) for i in labels],
        'chars': chars
    }
    os.mkdir('vis_cache')
    with open('vis_cache/images.npy', 'wb') as f:
        np.save(f, images)
    with open('vis_cache/km.npy', 'wb') as f:
        np.save(f, km)
    with open('vis_cache/meta.json', 'w') as f:
        json.dump(meta, f)


def coord_to_idx(a, b):
    return M_order * a + b

def idx_to_coord(idx):
    return (idx // M_order, idx % N_order)


if __name__ == '__main__':

    mask = labels < 20
    images_masked = images[mask]
    km_masked = km[mask]
    
    tsne = TSNE(n_components=2, verbose=0, perplexity=40, n_iter=300)
    tsne_res = tsne.fit_transform(km[mask])

    df = pd.DataFrame(tsne_res)
    df['label'] = pd.Series(labels[mask]).astype(str)
    df['index'] = df.index

    fig = px.scatter(df, x=0, y=1, color='label', labels={'label': 'Character'}, hover_data=['index'])
    fig.update_traces(
        hoverinfo="none",
        hovertemplate=None,
        marker=dict(size=10, opacity=0.8),
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)

    # update legend entries
    for i in range(len(fig.data)):
        fig.data[i].name = chars[int(fig.data[i].name)]

    app = Dash(__name__)


    app.layout = html.Div(
        className="container",
        children=[
            html.Div(children=[
                dcc.Graph(id="graph-2-dcc", figure=fig, clear_on_unhover=True),
                dcc.Tooltip(id="graph-tooltip-2", direction='bottom')
            ]),
            html.Div(id="div-nn", children=[], style={'padding-left': '70px', 'padding-right': '70px'})
        ],
        style={'font-family': 'arial', 'font-size': '11pt'}
    )

    @app.callback(
        Output("graph-tooltip-2", "show"),
        Output("graph-tooltip-2", "bbox"),
        Output("graph-tooltip-2", "children"),
        Output("graph-tooltip-2", "direction"),
        Input("graph-2-dcc", "hoverData"),
    )
    def display_hover(hoverData):
        
        if hoverData is None:
            return False, no_update, no_update, no_update

        hover_data = hoverData["points"][0]
        i = hover_data['customdata'][0]
        
        children = [
            html.Img(
                src=get_image_url(images_masked[i]),
                style={"width": "70px"},
            )
        ]

        direction = "bottom" if hover_data["y"] > 1.5 else "top"

        return True, hover_data["bbox"], children, direction


    @app.callback(
        Output('div-nn', 'children'),
        Input('graph-2-dcc', 'clickData')
    )
    def on_click(clickData):

        if clickData is None:
            return []

        click_data = clickData["points"][0]
        i = click_data['customdata'][0]
        km_masked[i]
        dists = ((km_masked - km_masked[i]) ** 2).sum(axis=1) / 1e7
        nns = dists.argsort()[1:] # nearest neighbours
        nns_class = nns[df['label'][nns] == df['label'][i]] # nn from same class
        ranks = np.array(range(len(nns))) + 1
        nn_elems = nns_to_elems(nns, dists, ranks)
        nn_class_elems = nns_to_elems(nns_class, dists, ranks[df['label'][nns] == df['label'][i]])
        div_style = {
            'white-space': 'nowrap', 
            'overflow-x': 'scroll', 
            'background-color': '#ededed',
            'padding': '20px',
            'margin-bottom': '20px'
        }
        return [
            html.Div(nn_elems, style=div_style),
            html.Div(nn_class_elems, style=div_style)
        ]


app.run_server(debug=True)








