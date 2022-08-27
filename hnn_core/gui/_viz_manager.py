"""HNN GUI visualization management tool"""

# Authors: Huzi Cheng <hzcheng15@icloud.com>

import copy
import io
from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from IPython.display import display
from ipywidgets import (Box, Button, Dropdown, FloatText, HBox, Label, Layout,
                        Output, Tab, VBox, link)

from hnn_core.gui._logging import logger
from hnn_core.viz import plot_dipole


class _TabWithNoneIndex(Tab):
    def _reset_selected_index(self):
        """Disable the default tab behavior that reset selected_index to 0."""
        # if there are no tabs, then none should be selected
        num_children = len(self.children)
        if num_children == 0:
            self.selected_index = None

        # if there are tabs, but none is selected, select the first one
        elif self.selected_index is None:
            self.selected_index = None

        # if there are tabs and a selection, but the selection is no longer
        # valid, select the last tab.
        elif num_children < self.selected_index:
            self.selected_index = num_children - 1


_dpi = 96
_fig_placeholder = 'Run simulation to add figures here.'

_plot_types = [
    'current dipole',
    'dipole@L2',
    'dipole@L5',
    'dipole@agg',
    'input histogram',
    'spikes',
    'PSD',
    'spectogram',
    'network',
    'layer-specific dipole',
]

_spectrogram_color_maps = [
    "viridis",
    "plasma",
    "inferno",
    "magma",
    "cividis",
]


fig_templates = {
    "2row x 1col (1:3)": {
        "kwargs":
        "gridspec_kw={\"height_ratios\":[1,3]}",
        "mosaic": "00\n11",
    },
    "2row x 1col (1:1)": {
        "kwargs":
        "gridspec_kw={\"height_ratios\":[1,1]}",
        "mosaic": "00\n11",
    },
    "1row x 2col (1:1)": {
        "kwargs":
        "gridspec_kw={\"height_ratios\":[1,1]}",
        "mosaic": "01\n01",
    },
    "single figure": {
        "kwargs": "",
        "mosaic": "00\n00",
    },
    "2row x 2col (1:1)": {
        "kwargs":
        "gridspec_kw={\"height_ratios\":[1,1]}",
        "mosaic": "01\n23",
    },
}


def _idx2figname(idx):
    return f"Figure {idx}"


def _figname2idx(fname):
    return int(fname.split(" ")[-1])


def _update_ax(fig, ax, single_simulation, plot_type, plot_config):
    """Refresh plots with simulation_data.

    Parameters
    ----------
    fig: Figure
        A matplotlib.figure.Figure object.
    ax: Axes
        matplotlib.axes.Axes
    single_simulation: dict
        A single simulation
    plot_type: str
        Type of subplots
    plot_config: dict
        A dict specifies the preprocessing and style of plots.
    """
    # Make sure that visualization does not change the original data
    dpls_copied = copy.deepcopy(single_simulation['dpls'])
    net_copied = copy.deepcopy(single_simulation['net'])

    for dpl in dpls_copied:
        dpl.smooth(plot_config['dipole_smooth']).scale(
            plot_config['dipole_scaling'])

    if net_copied is None:
        print("No network data")
        return

    if plot_type == 'spikes':
        if net_copied.cell_response:
            net_copied.cell_response.plot_spikes_raster(ax=ax, show=False)

    elif plot_type == 'input histogram':
        # BUG: got error here, need a better way to handle exception
        if net_copied.cell_response:
            net_copied.cell_response.plot_spikes_hist(ax=ax, show=False)

    elif plot_type == 'PSD':
        if len(dpls_copied) > 0:
            dpls_copied[0].plot_psd(fmin=0, fmax=50, ax=ax, show=False)

    elif plot_type == 'spectogram':
        if len(dpls_copied) > 0:
            min_f = 10.0
            max_f = plot_config['max_spectral_frequency']
            step_f = 1.0
            freqs = np.arange(min_f, max_f, step_f)
            n_cycles = freqs / 8.
            dpls_copied[0].plot_tfr_morlet(
                freqs,
                n_cycles=n_cycles,
                colormap=plot_config['spectrogram_cm'],
                ax=ax, show=False)

    elif 'dipole@' in plot_type:
        if len(dpls_copied) > 0:
            plot_dipole(dpls_copied,
                        ax=ax,
                        layer=plot_type.split("@")[1],
                        average=True,
                        show=False)

    elif plot_type == 'current dipole':
        if len(dpls_copied) > 0:
            plot_dipole(dpls_copied, ax=ax, average=True, show=False)

    elif plot_type == 'layer-specific dipole':
        if len(dpls_copied) > 0:
            ax.remove()
            layers = ["L2", "L5", "agg"]
            gridspec = ax.get_subplotspec()
            gs01 = gridspec.subgridspec(3, 1)
            ax_l2 = fig.add_subplot(gs01[0])
            ax_l5 = fig.add_subplot(gs01[1])
            ax_lagg = fig.add_subplot(gs01[2])
            axes = [ax_l2, ax_l5, ax_lagg]
            plot_dipole(dpls_copied, ax=axes,
                        layer=layers, average=True)
        else:
            print("No dipole data")

    elif plot_type == 'network':
        # raise NotImplementedError("No 3D plot for single ax now")
        if net_copied:
            _fig = plt.figure()
            _ax = _fig.add_subplot(111, projection='3d')
            net_copied.plot_cells(ax=_ax)

            io_buf = io.BytesIO()
            _fig.savefig(io_buf, format='raw')
            io_buf.seek(0)
            img_arr = np.reshape(np.frombuffer(io_buf.getvalue(),
                                               dtype=np.uint8),
                                 newshape=(int(_fig.bbox.bounds[3]),
                                           int(_fig.bbox.bounds[2]), -1))
            io_buf.close()
            ax.imshow(img_arr)


def _plot_on_axes(b, widgets_simulation, widgets_plot_type,
                  spectrogram_colormap_selection, dipole_smooth,
                  max_spectral_frequency, dipole_scaling, data, fig, ax):
    single_simulation = data['simulations'][widgets_simulation.value]

    plot_config = {
        "max_spectral_frequency": max_spectral_frequency.value,
        "dipole_scaling": dipole_scaling.value,
        "dipole_smooth": dipole_smooth.value,
        "spectrogram_cm": spectrogram_colormap_selection.value
    }
    plot_type = widgets_plot_type.value
    ax.clear()
    ax.set_facecolor('w')
    _update_ax(fig, ax, single_simulation, plot_type, plot_config)


def _clear_axis(b, ax):
    ax.clear()
    ax.set_facecolor('w')


def _get_ax_control(widgets, data, fig_idx, fig, ax):
    analysis_style = {'description_width': '200px'}
    layout = Layout(width="98%")
    simulation_names = tuple(data['simulations'].keys())
    if len(simulation_names) == 0:
        simulation_names = ["None", ]

    simulation_selection = Dropdown(
        options=simulation_names,
        value=simulation_names[-1], description='Simulation:',
        disabled=False, layout=layout,
        style=analysis_style,
    )

    plot_type_selection = Dropdown(
        options=_plot_types,
        value=_plot_types[0], description='Type:',
        disabled=False, layout=layout,
        style=analysis_style,
    )

    spectrogram_colormap_selection = Dropdown(
        description='Spectrogram Colormap:',
        options=[(cm, cm) for cm in _spectrogram_color_maps],
        value=_spectrogram_color_maps[0], layout=layout,
        style=analysis_style,
    )
    dipole_smooth = FloatText(value=30,
                              description='Dipole Smooth Window (ms):',
                              disabled=False, layout=layout,
                              style=analysis_style)
    dipole_scaling = FloatText(
        value=3000,
        description='Dipole Scaling:',
        disabled=False, layout=layout,
        style=analysis_style)

    max_spectral_frequency = FloatText(
        value=100,
        description='Max Spectral Frequency (Hz):',
        disabled=False, layout=layout,
        style=analysis_style)

    clear_button = Button(description='Clear axis')
    clear_button.on_click(partial(_clear_axis, ax=ax))

    plot_button = Button(description='Plot')
    plot_button.on_click(
        partial(
            _plot_on_axes,
            widgets_simulation=simulation_selection,
            widgets_plot_type=plot_type_selection,
            spectrogram_colormap_selection=spectrogram_colormap_selection,
            dipole_smooth=dipole_smooth,
            max_spectral_frequency=max_spectral_frequency,
            dipole_scaling=dipole_scaling,
            data=data,
            fig=fig,
            ax=ax,
        ))

    vbox = VBox([
        simulation_selection, plot_type_selection, dipole_smooth,
        dipole_scaling, max_spectral_frequency, spectrogram_colormap_selection,
        HBox(
            [plot_button, clear_button],
            layout=Layout(justify_content='space-between', ),
        )
    ], layout=Layout(width="98%"))

    return vbox


def _close_figure(b, widgets, data, fig_idx):
    fig_related_widgets = [widgets['figs_tabs'], widgets['axes_config_tabs']]
    for w_idx, tab in enumerate(fig_related_widgets):
        tab_children = list(tab.children)
        titles = [tab.get_title(idx) for idx in range(len(tab.children))]
        tab_idx = titles.index(_idx2figname(fig_idx))
        print(f"Del fig_idx={fig_idx}, fig_idx={fig_idx}")
        del tab_children[tab_idx], titles[tab_idx]

        tab.children = tuple(tab_children)
        [tab.set_title(idx, title) for idx, title in enumerate(titles)]

        if w_idx == 0:
            plt.close(data['figs'][fig_idx])
            del data['figs'][fig_idx]
            n_tabs = len(tab.children)
            for idx in range(n_tabs):
                _fig_idx = _figname2idx(tab.get_title(idx))
                assert _fig_idx in data['figs'].keys()

                tab.children[idx].clear_output()
                with tab.children[idx]:
                    display(data['figs'][_fig_idx].canvas)

        if n_tabs == 0:
            widgets['figs_output'].clear_output()
            with widgets['figs_output']:
                display(Label(_fig_placeholder))


def _add_axes_controls(widgets, data, fig, axd):
    fig_idx = data['fig_idx']['idx']

    controls = Tab()
    children = [
        _get_ax_control(widgets, data, fig_idx=fig_idx, fig=fig, ax=ax)
        for ax_key, ax in axd.items()
    ]
    controls.children = children
    for i in range(len(children)):
        controls.set_title(i, f'ax{i}')

    close_fig_button = Button(description=f'Close {_idx2figname(fig_idx)}',
                              button_style='danger', icon='close',
                              layout=Layout(width="98%"))
    close_fig_button.on_click(
        partial(_close_figure, widgets=widgets, data=data, fig_idx=fig_idx))

    _prev_titles = copy.deepcopy(widgets['axes_config_tabs'].titles)
    widgets['axes_config_tabs'].children = widgets[
        'axes_config_tabs'].children + (VBox([close_fig_button, controls]), )
    widgets['axes_config_tabs'].titles = _prev_titles + (
        _idx2figname(fig_idx), )


def _add_figure(b, widgets, data, scale=0.95):
    template_name = widgets['templates_dropdown'].value
    fig_idx = data['fig_idx']['idx']
    viz_output_layout = data['visualization_output']
    fig_outputs = Output()
    n_tabs = len(widgets['figs_tabs'].children)
    if n_tabs == 0:
        widgets['figs_output'].clear_output()
        with widgets['figs_output']:
            display(widgets['figs_tabs'])

    _prev_titles = copy.deepcopy(widgets['figs_tabs'].titles)
    widgets['figs_tabs'].selected_index = None
    widgets['figs_tabs'].children = widgets['figs_tabs'].children + (
        fig_outputs, )

    widgets['figs_tabs'].titles = _prev_titles + (_idx2figname(fig_idx), )

    with fig_outputs:
        dpi = _dpi
        figsize = (scale * ((int(viz_output_layout.width[:-2]) - 10) / dpi),
                   scale * ((int(viz_output_layout.height[:-2]) - 10) / dpi))
        mosaic = fig_templates[template_name]['mosaic']
        kwargs = eval(f"dict({fig_templates[template_name]['kwargs']})")
        fig, axd = plt.subplot_mosaic(mosaic, figsize=figsize, dpi=dpi,
                                      **kwargs)
        fig.tight_layout()
        fig.canvas.header_visible = False
        fig.canvas.footer_visible = False
        plt.show()

    _add_axes_controls(widgets, data, fig=fig, axd=axd)

    data['figs'][fig_idx] = fig
    data['fig_idx']['idx'] += 1
    widgets['figs_tabs'].selected_index = n_tabs


class _VizManager:
    """GUI visualization panel manager class.

    Parameters
    ----------
    gui_data : dict
        A dict containing all simulation data
    viz_layout : dict
        A dict about visualization layout specs
    """
    def __init__(self, gui_data, viz_layout):
        plt.close("all")
        self.viz_layout = viz_layout

        self.axes_config_output = Output()
        self.figs_output = Output(
            layout=self.viz_layout['visualization_window'])

        # widgets
        self.axes_config_tabs = _TabWithNoneIndex()
        self.figs_tabs = _TabWithNoneIndex()
        self.axes_config_tabs.selected_index = None
        self.figs_tabs.selected_index = None
        link(
            (self.axes_config_tabs, 'selected_index'),
            (self.figs_tabs, 'selected_index'),
        )

        template_names = list(fig_templates.keys())
        self.templates_dropdown = Dropdown(
            description='Layout template:',
            options=template_names,
            value=template_names[0],
            style={'description_width': 'initial'},
            layout=Layout(width="98%"))
        self.make_fig_button = Button(
            description='Make figure',
            button_style="primary",
            style={'button_color': self.viz_layout['theme_color']},
            layout=self.viz_layout['btn'])
        self.make_fig_button.on_click(self.add_figure)

        # data
        self.fig_idx = {"idx": 1}
        self.figs = {}
        self.gui_data = gui_data

    @property
    def widgets(self):
        return {
            "figs_output": self.figs_output,
            "axes_config_tabs": self.axes_config_tabs,
            "figs_tabs": self.figs_tabs,
            "templates_dropdown": self.templates_dropdown
        }

    @property
    def data(self):
        return {
            "simulations": self.gui_data["simulation_data"],
            "fig_idx": self.fig_idx,
            "visualization_output": self.viz_layout['visualization_output'],
            "figs": self.figs
        }

    def compose(self):
        with self.axes_config_output:
            display(self.axes_config_tabs)
        with self.figs_output:
            display(Label(_fig_placeholder))

        config_panel = VBox([
            Box(
                [
                    self.templates_dropdown,
                    self.make_fig_button,
                ],
                layout=Layout(
                    display='flex',
                    flex_flow='column',
                    align_items='stretch',
                ),
            ),
            Label("Figure config:"),
            self.axes_config_output,
        ])
        return config_panel, self.figs_output

    def add_figure(self, b=None):
        logger.debug("add figure")
        _add_figure(b, self.widgets, self.data, scale=0.97)