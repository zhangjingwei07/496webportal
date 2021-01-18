def fig_write_return(fig, ret_type, save):
    if save:
        fig.write_image(save)
    if ret_type == "fig":
        return fig
    else:
        return fig.to_json()
