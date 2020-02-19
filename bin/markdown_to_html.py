#!/usr/bin/env python
from __future__ import print_function
import argparse
import markdown
import os
import sys

def convert_markdown(in_fn):
    input_md = open(in_fn, mode="r", encoding="utf-8").read()
    html = markdown.markdown(
        "[TOC]\n" + input_md,
        extensions = [
            'pymdownx.extra',
            'pymdownx.b64',
            'pymdownx.highlight',
            'pymdownx.emoji',
            'pymdownx.tilde',
            'toc'
        ],
        extension_configs = {
            'pymdownx.b64': {
                'base_path': os.path.dirname(in_fn)
            },
            'pymdownx.highlight': {
                'noclasses': True
            },
            'toc': {
                'title': 'Table of Contents'
            }
        }
    )
    return html

def wrap_html(contents):
    header = """<!DOCTYPE html><html>
    <head>
        <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
        <style>
            body {
              font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,"Noto Sans",sans-serif,"Apple Color Emoji","Segoe UI Emoji","Segoe UI Symbol","Noto Color Emoji";
              padding: 3em;
              margin-right: 350px;
              max-width: 100%;
            }
            .toc {
              position: fixed;
              right: 20px;
              width: 300px;
              padding-top: 20px;
              overflow: scroll;
              height: calc(100% - 3em - 20px);
            }
            .toctitle {
              font-size: 1.8em;
              font-weight: bold;
            }
            .toc > ul {
              padding: 0;
              margin: 1rem 0;
              list-style-type: none;
            }
            .toc > ul ul { padding-left: 20px; }
            .toc > ul > li > a { display: none; }
            img { max-width: 800px; }
            pre {
              padding: 0.6em 1em;
            }
            h2 {

            }
        </style>
    </head>
    <body>
    <div class="container">
    """
    footer = """
    </div>
    </body>
    </html>
    """
    return header + contents + footer


def parse_args(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('mdfile', type=argparse.FileType('r'), nargs='?',
                        help='File to convert. Defaults to stdin.')
    parser.add_argument('-o', '--out', type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='Output file name. Defaults to stdout.')
    return parser.parse_args(args)

def main(args=None):
    args = parse_args(args)
    converted_md = convert_markdown(args.mdfile.name)
    html = wrap_html(converted_md)
    args.out.write(html)

if __name__ == '__main__':
    sys.exit(main())
