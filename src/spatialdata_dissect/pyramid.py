"""Access to image pyramid levels inside a SpatialData object."""


def pyramid_levels(sdata, image_key="morphology_focus"):
    """Return image pyramid levels as DataArrays, sorted smallest-first."""
    image = sdata.images[image_key]
    if hasattr(image, "shape") and hasattr(image, "dims"):
        return [image]

    levels = []
    for child in getattr(image, "children", {}).values():
        dataset = child.to_dataset()
        if not dataset.data_vars:
            continue
        name = "image" if "image" in dataset.data_vars else next(iter(dataset.data_vars))
        levels.append(dataset[name])

    if not levels:
        raise ValueError(f"No image levels found for {image_key!r}")

    return sorted(levels, key=lambda a: int(a.sizes["y"]) * int(a.sizes["x"]))
