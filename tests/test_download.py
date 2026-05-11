"""Offline tests for the bulk-download retry classifier."""

from __future__ import annotations

import tarfile

import pytest
import requests

from tcga_pull.download import _is_retriable_http_error


def _http_error(status: int) -> requests.HTTPError:
    resp = requests.Response()
    resp.status_code = status
    return requests.HTTPError(f"HTTP {status}", response=resp)


@pytest.mark.parametrize("status", [500, 502, 503, 504])
def test_retriable_5xx(status: int):
    assert _is_retriable_http_error(_http_error(status))


@pytest.mark.parametrize("status", [400, 401, 403, 404, 422])
def test_not_retriable_4xx(status: int):
    assert not _is_retriable_http_error(_http_error(status))


def test_retriable_connection_error():
    assert _is_retriable_http_error(requests.ConnectionError("network down"))


def test_retriable_timeout():
    assert _is_retriable_http_error(requests.Timeout("read timed out"))


def test_retriable_tar_error():
    assert _is_retriable_http_error(tarfile.TarError("bad tar"))


def test_not_retriable_arbitrary_exception():
    assert not _is_retriable_http_error(ValueError("nope"))
    assert not _is_retriable_http_error(KeyError("missing"))
