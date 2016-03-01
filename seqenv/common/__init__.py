# Built-in modules #
import tempfile, datetime, dateutil, hashlib

# One liners #
flatter = lambda x: [item for sublist in x for item in sublist]

################################################################################
def new_temp_path(**kwargs):
    """A new temporary path"""
    handle = tempfile.NamedTemporaryFile(**kwargs)
    path = handle.name
    handle.close()
    return path

################################################################################
class GenWithLength(object):
    """A generator with a length attribute"""
    def __init__(self, gen, length): self.gen, self.length = gen, length
    def __iter__(self): return self.gen
    def __len__(self): return self.length

################################################################################
def pretty_now():
    """Returns some thing like '2014-07-24 11:12:45 CEST+0200'"""
    now = datetime.datetime.now(dateutil.tz.tzlocal())
    return now.strftime("%Y-%m-%d %H:%M:%S %Z%z")

################################################################################
def download_from_url(source, destination, progress=False):
    """Download a file from an URL and place it somewhere. Like wget.
    Uses requests and tqdm to display progress"""
    from tqdm import tqdm
    import requests
    response = requests.get(source, stream=True)
    with open(destination, "wb") as handle:
        if progress:
            for data in tqdm(response.iter_content()): handle.write(data)
        else:
            for data in response.iter_content(): handle.write(data)

################################################################################
def md5sum(file_path, blocksize=65536):
    """Compute the md5 of a file. Pretty fast."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(blocksize), ""):
            md5.update(block)
    return md5.hexdigest()