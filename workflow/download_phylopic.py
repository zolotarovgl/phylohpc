#!/usr/bin/env python3
"""Search PhyloPic by taxon name and download a 192x192 thumbnail."""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.parse
import urllib.request
import urllib.error
from pathlib import Path
from typing import Any


API_BASE = "https://api.phylopic.org"
SITE_BASE = "https://www.phylopic.org"
USER_AGENT = "phylohpc-phylopic-downloader/1.0"


def fetch_json(url: str) -> dict[str, Any]:
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urllib.request.urlopen(req) as resp:
            return json.load(resp)
    except urllib.error.HTTPError as exc:
        body = exc.read()
        try:
            return json.loads(body.decode("utf-8"))
        except Exception as decode_exc:  # pragma: no cover - best-effort error detail
            raise RuntimeError(f"HTTP {exc.code} for {url}") from decode_exc


def fetch_bytes(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as resp:
        return resp.read()


def fetch_text(url: str) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as resp:
        return resp.read().decode("utf-8")


def absolutize_href(href: str) -> str:
    if href.startswith("http://") or href.startswith("https://"):
        return href
    return urllib.parse.urljoin(API_BASE, href)


def quote_query(name: str) -> str:
    return urllib.parse.quote(name.strip(), safe="")


def list_matching_nodes(name: str) -> list[dict[str, str]]:
    search_url = f"{SITE_BASE}/search?q={quote_query(name)}"
    html = fetch_text(search_url)
    raw_matches = re.findall(r'"/nodes/([0-9a-f-]+)\?build=(\d+)","([^"]+)"', html, flags=re.I)
    out: list[dict[str, str]] = []
    seen: set[str] = set()
    for node_uuid, build, title in raw_matches:
        href = f"{API_BASE}/nodes/{node_uuid}?build={build}"
        if href in seen:
            continue
        seen.add(href)
        out.append({"title": title, "href": href})
    return out


def choose_match(matches: list[dict[str, str]], name: str, index: int | None) -> tuple[dict[str, str], bool]:
    if not matches:
        raise RuntimeError(f"No PhyloPic matches found for {name!r}.")
    if index is not None:
        if index < 0 or index >= len(matches):
            raise RuntimeError(f"Requested result index {index} is out of range for {len(matches)} matches.")
        chosen = matches[index]
        return chosen, chosen["title"].casefold() == name.casefold()
    wanted = name.casefold()
    for match in matches:
        if match["title"].casefold() == wanted:
            return match, True
    return matches[0], False


def get_primary_image_url(node_url: str) -> str:
    node = fetch_json(node_url)
    image_href = node.get("_links", {}).get("primaryImage", {}).get("href")
    if not image_href:
        raise RuntimeError(f"No primary image is available for node {node_url}.")
    return absolutize_href(image_href)


def pick_thumbnail_url(image_url: str, size: str = "192x192") -> tuple[str, dict[str, Any]]:
    image = fetch_json(image_url)
    thumbs = image.get("_links", {}).get("thumbnailFiles", [])
    if not thumbs:
        raise RuntimeError(f"No thumbnail files are available for image {image_url}.")
    for thumb in thumbs:
        if thumb.get("sizes") == size and thumb.get("href"):
            return str(thumb["href"]), image
    for thumb in thumbs:
        if thumb.get("href"):
            return str(thumb["href"]), image
    raise RuntimeError(f"Could not find a usable thumbnail URL for image {image_url}.")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Search PhyloPic for a taxon name and download a 192x192 thumbnail PNG."
    )
    p.add_argument("name", help="Species or taxon name to search for, e.g. 'Homo sapiens'.")
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output PNG path. Defaults to ./<species_name>.png.",
    )
    p.add_argument(
        "--result-index",
        type=int,
        default=None,
        help="0-based search-result index to use. By default the first exact title match is preferred, otherwise result 0.",
    )
    p.add_argument(
        "--list",
        action="store_true",
        help="List matching PhyloPic nodes and exit without downloading.",
    )
    return p


def default_output_path(name: str) -> Path:
    slug = "_".join(name.strip().split())
    return Path(f"{slug}.png")


def main() -> int:
    args = build_parser().parse_args()
    try:
        matches = list_matching_nodes(args.name)
        if args.list:
            if not matches:
                print(f"No PhyloPic matches found for {args.name!r}.", file=sys.stderr)
                return 1
            for i, match in enumerate(matches):
                print(f"[{i}] {match['title']}\t{match['href']}")
            return 0

        match, exact = choose_match(matches, args.name, args.result_index)
        image_url = get_primary_image_url(match["href"])
        thumb_url, image_meta = pick_thumbnail_url(image_url, "192x192")
        output = args.output or default_output_path(args.name)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(fetch_bytes(thumb_url))

        image_title = image_meta.get("_links", {}).get("self", {}).get("title", "")
        attribution = image_meta.get("attribution", "")
        license_url = image_meta.get("_links", {}).get("license", {}).get("href", "")

        if not exact and args.result_index is None:
            print(
                f"Warning: no exact title match for {args.name!r}; using first result {match['title']!r}.",
                file=sys.stderr,
            )
        print(f"Saved {output}")
        print(f"Node: {match['title']}")
        if image_title:
            print(f"Image: {image_title}")
        if attribution:
            print(f"Attribution: {attribution}")
        if license_url:
            print(f"License: {license_url}")
        print(f"Thumbnail: {thumb_url}")
        return 0
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
