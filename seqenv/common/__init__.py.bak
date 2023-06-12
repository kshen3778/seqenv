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
def md5sum(file_path, blocksize=65536):
    """Compute the md5 of a file. Pretty fast."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for block in iter(lambda: f.read(blocksize), ""):
            md5.update(block)
    return md5.hexdigest()

################################################################################
def download_from_url(source, destination, progress=False, uncompress=True):
    """Download a file from an URL and place it somewhere. Like wget.
    Uses requests and tqdm to display progress if you want.
    By default it will uncompress files."""
    # Modules #
    from tqdm import tqdm
    import requests
    from autopaths import FilePath
    # Check destination exists #
    destination = FilePath(destination)
    destination.directory.create_if_not_exists()
    # Over HTTP #
    response = requests.get(source, stream=True)
    total_size = int(response.headers.get('content-length'))
    block_size = total_size/1024
    # Do it #
    with open(destination, "wb") as handle:
        if progress:
            for data in tqdm(response.iter_content(chunk_size=block_size), total=1024): handle.write(data)
        else:
            for data in response.iter_content(chunk_size=block_size): handle.write(data)
    # Uncompress #
    if uncompress:
        with open(destination) as f: header = f.read(4)
        if header == "PK\x03\x04": unzip(destination, inplace=True)
        # Add other compression formats here
    # Return #
    return destination

################################################################################
def unzip(source, destination=None, inplace=False, single=True):
    """Unzip a standard zip file. Can specify the destination of the
    uncompressed file, or just set inplace=True to delete the original."""
    # Load #
    import zipfile, tempfile, shutil
    # Check #
    assert zipfile.is_zipfile(source)
    # Load #
    z = zipfile.ZipFile(source)
    if single or inplace: assert len(z.infolist()) == 1
    # Single file #
    if single:
        member = z.infolist()[0]
        tmpdir = tempfile.mkdtemp() + '/'
        z.extract(member, tmpdir)
        z.close()
        if inplace: shutil.move(tmpdir + member.filename, source)
        else:       shutil.move(tmpdir + member.filename, destination)
    # Multifile - no security, dangerous #
    if not single:
        z.extractall()
