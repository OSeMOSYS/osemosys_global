"""Functions to set storage levels."""

def set_storage_level_start(storage_set, region_name):

    storage_level_start = storage_set.copy().rename(columns = {'VALUE' : 'STORAGE'})
    storage_level_start.insert(loc=0, column = 'REGION', value = region_name)
    storage_level_start.insert(loc=2, column = 'VALUE', value = 0)

    return storage_level_start