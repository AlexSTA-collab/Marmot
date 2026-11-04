"""Utilities for comparing traces on vertically extruded inclined layers."""

import numpy as np


def interpolate_top_trace_along_normal(
    points_bottom,
    points_top,
    values_top,
    geometry_height,
    angle_deg,
    tol=1.0e-6,
):
    """Interpolate the upper trace at normal correspondents of lower points.

    The mesh generator offsets the upper face by ``(0, 0, H)``.  For a face
    inclined by ``angle_deg``, the normal separation is ``H cos(angle)`` and
    normal-corresponding points have a nonzero in-plane offset.  ``values_top``
    may be a vector, tensor, or any array whose first axis indexes points.

    Returns the valid lower-point mask, the corresponding upper points, the
    interpolated upper values, and the normal offset vector.
    """
    points_bottom = np.asarray(points_bottom, dtype=float)
    points_top = np.asarray(points_top, dtype=float)
    values_top = np.asarray(values_top)

    theta = np.radians(float(angle_deg))
    normal = np.array([np.sin(theta), 0.0, np.cos(theta)])
    normal_offset = float(geometry_height) * np.cos(theta) * normal
    query_points = points_bottom + normal_offset

    # The extracted nodal and cell-center traces are tensor-product XY grids.
    # Integer keys avoid fragile direct floating-point equality checks.
    x_keys = np.rint(points_top[:, 0] / tol).astype(np.int64)
    y_keys = np.rint(points_top[:, 1] / tol).astype(np.int64)
    unique_x_keys = np.unique(x_keys)
    unique_y_keys = np.unique(y_keys)
    nx, ny = len(unique_x_keys), len(unique_y_keys)
    if nx * ny != len(points_top):
        raise ValueError(
            "Upper trace is not a complete tensor-product XY grid: "
            f"{len(points_top)} points, expected {nx} x {ny} = {nx * ny}."
        )

    x_index = {key: i for i, key in enumerate(unique_x_keys)}
    y_index = {key: i for i, key in enumerate(unique_y_keys)}
    value_shape = values_top.shape[1:]
    grid = np.empty((ny, nx) + value_shape, dtype=values_top.dtype)
    for point_index, (x_key, y_key) in enumerate(zip(x_keys, y_keys)):
        grid[y_index[y_key], x_index[x_key]] = values_top[point_index]

    xs = unique_x_keys.astype(float) * tol
    ys = unique_y_keys.astype(float) * tol
    qx = query_points[:, 0]
    qy = query_points[:, 1]
    valid = (
        (qx >= xs[0] - tol)
        & (qx <= xs[-1] + tol)
        & (qy >= ys[0] - tol)
        & (qy <= ys[-1] + tol)
    )

    qx_valid = np.clip(qx[valid], xs[0], xs[-1])
    qy_valid = np.clip(qy[valid], ys[0], ys[-1])
    ix = np.clip(np.searchsorted(xs, qx_valid, side="right") - 1, 0, nx - 2)
    iy = np.clip(np.searchsorted(ys, qy_valid, side="right") - 1, 0, ny - 2)

    wx = (qx_valid - xs[ix]) / (xs[ix + 1] - xs[ix])
    wy = (qy_valid - ys[iy]) / (ys[iy + 1] - ys[iy])
    weight_shape = (len(wx),) + (1,) * len(value_shape)
    wx = wx.reshape(weight_shape)
    wy = wy.reshape(weight_shape)

    v00 = grid[iy, ix]
    v10 = grid[iy, ix + 1]
    v01 = grid[iy + 1, ix]
    v11 = grid[iy + 1, ix + 1]
    interpolated = (
        (1.0 - wx) * (1.0 - wy) * v00
        + wx * (1.0 - wy) * v10
        + (1.0 - wx) * wy * v01
        + wx * wy * v11
    )

    return valid, query_points[valid], interpolated, normal_offset
