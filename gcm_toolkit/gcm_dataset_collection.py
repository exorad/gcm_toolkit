""" GCM dataset collection class to deal with GCM data """
from collections import UserDict
from .core import writer as wrt


class GCMDatasetCollection(UserDict):
    """
    This class represents a collection of 3D GCM Datasets.
    A GCMDatasetCollection is a dictionary with in which GCM models are loaded
    with a tag.
    """

    def get_models(self, tag=None, always_dict=False):
        """
        Function return all GCMs in memory. If a tag is given, only return this
        one.

        Parameters
        ----------
        tag : str
            Name of the model that should be returned.
        always_dict: bool
            Force result to be a dictionary (if tag is None)

        Returns
        -------
        selected_models : GCMDatasetCollection or xarray Dataset
            All models in self._models, or only the one with the right tag.
            Will definetly be GCMDatasetCollection if always_dict=True
        """
        if isinstance(tag, str):
            return self.get(tag)

        # If no tag is given, return all models
        if tag is None:
            if len(self) > 1 or always_dict:
                return self
            return list(self.values())[0]

        # If the tag is not a string, raise an error
        return wrt.write_status("ERROR", "The given tag is not a string.")

    def get_one_model(self, tag=None, raise_error=True):
        """
        Helper Function that raises an error or returns None,
        if more than one model is selected.

        Parameters
        ----------
        tag: str, optional
            Name of the model that should be returned
        raise_error: bool, optional
            If true, function will raise error,
            else will return None if not one model is selected

        Returns
        -------
        ds: xarray Dataset
            Selected model
        """

        # select the appropriate dataset
        dsi = self.get_models(tag=tag, always_dict=False)
        # Raise error, if key not in collection and raise error is specified
        if dsi is None and raise_error:
            raise KeyError("No dataset for given key available.")
        # if a collection is given (because multiple datasets are available, and
        # the tag is not provided), avoid ambiguity by raising an error
        if isinstance(dsi, GCMDatasetCollection):
            wrt.write_status("ERROR", "Ambiguous task. Please provide a tag.")

        return dsi
