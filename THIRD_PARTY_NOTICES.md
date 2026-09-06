# Third Party Notices

## Vendored mapper seeding code

The vendored mapper crate under `crates/metam` contains seed-generation code
that substantially adapts Strobealign-style syncmer and randstrobe seeding.
The source file `crates/metam/src/seeds.rs` includes a source-level
attribution comment.

Strobealign is distributed under the MIT license. The MIT license text follows:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Other dependencies

All other dependencies are obtained from crates.io at build time under their
own licences (MIT, Apache-2.0, BSD, ISC, Zlib, Unicode-3.0, and MPL-2.0 for
`option-ext`). Their exact versions are recorded in `Cargo.lock`, and each
crate's licence text is included in its published source. The bundled htslib
sources compiled by `rust-htslib` are distributed under the MIT licence.
