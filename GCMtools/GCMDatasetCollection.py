import GCMtools.core.writer as wrt
from collections import UserDict


class GCMDatasetCollection(UserDict):
    """
    This class represents a collection of 3D GCM Datasets.
    A GCMDatasetCollection is a dictionary with in which GCM models are loaded
    with a tag.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_models(self, tag=None):
        """
        Function return all GCMs in memory. If a tag is given, only return this
        one.

        Parameters
        ----------
        tag : str
            Name of the model that should be returned.

        Returns
        -------
        selected_models : GCMDatasetCollection or xarray Dataset
            All models in self._models, or only the one with the right tag.
        """
        # If no tag is given, return all models
        if tag is None:
            return self._prepare_models()
        # If the tag is valid, return the corresponding model
        elif isinstance(tag, str):
            return self[tag]
        # If the tag is not a string, raise an error
        else:
            wrt.write_status('ERROR', 'The given tag is not a string.')

    def get_one_model(self, tag=None):
        """
        Helper Function that raises an error, if more than one model is selected.

        Parameters
        ----------
        tag: str, optional
            Name of the model that should be returned
        Returns
        -------
        ds: xarray Dataset
            Selected model
        """

        # select the appropriate dataset
        ds = self.get_models(tag=tag)
        # if a collection is given (because multiple datasets are available, and
        # the tag is not provided), avoid ambiguity by raising an error
        if isinstance(ds, GCMDatasetCollection) and len(ds) > 1 and tag is None:
            wrt.write_status('ERROR', 'Ambiguous task. Please provide a tag.')

    def _prepare_models(self):
        """
        returns dataset if only one model is loaded, otherwise collection of datasets.
        """
        if len(self) > 1:
            return self
        else:
            return list(self.values())[0]
