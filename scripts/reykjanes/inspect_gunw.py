"""Open one NISAR GUNW granule remotely (HTTP range requests) and list its layout."""
import sys
from pathlib import Path

import fsspec
import h5py

ROOT = Path(__file__).resolve().parents[2]


def token():
    from _secrets import read_env
    return read_env()["EARTHDATA_API_KEY"]


def presigned(url):
    """ASF answers every request with a redirect to a signed CloudFront URL (valid ~1 h).
    Resolving it once avoids paying the redirect on every range request."""
    import requests
    assert url.startswith("https://nisar.asf.earthdatacloud.nasa.gov/")
    r = requests.get(url, headers={"Authorization": f"Bearer {token()}", "Range": "bytes=0-0"},
                     allow_redirects=False, timeout=60)
    loc = r.headers.get("Location")
    if r.status_code in (301, 302, 303, 307, 308) and loc and loc.startswith("https://"):
        return loc, {}
    return url, {"Authorization": f"Bearer {token()}"}


def open_remote(url, block_mb=4):
    href, headers = presigned(url)
    import aiohttp
    fs = fsspec.filesystem("https", headers=headers, block_size=block_mb * 1024 * 1024,
                           client_kwargs={"timeout": aiohttp.ClientTimeout(total=None, connect=60, sock_read=180)})
    return h5py.File(fs.open(href, "rb", cache_type="blockcache"), "r")


if __name__ == "__main__":
    url = sys.argv[1]
    prefix = sys.argv[2] if len(sys.argv) > 2 else "/"
    with open_remote(url) as f:
        print("size GB", f.id.get_filesize() / 1e9)

        def show(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(name, obj.shape, obj.dtype, obj.chunks, obj.compression,
                      f"{obj.id.get_storage_size() / 1e6:.1f}MB")
        f[prefix].visititems(show)
