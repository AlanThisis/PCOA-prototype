from __future__ import annotations

import pandas as pd

from diversity import compute_pcoa


def test_compute_pcoa_returns_expected_shapes_for_braycurtis() -> None:
    feature_table = pd.DataFrame(
        {
            "feat1": [10, 0, 4, 2],
            "feat2": [0, 7, 1, 3],
            "feat3": [5, 1, 0, 8],
        },
        index=["sample1", "sample2", "sample3", "sample4"],
    )

    dm, ord_res = compute_pcoa(feature_table, metric="braycurtis")

    assert dm.shape == (4, 4)
    assert list(dm.ids) == ["sample1", "sample2", "sample3", "sample4"]
    assert ord_res.samples.shape[0] == 4
