---
title: "Site relaunch"
description: "A short note about why I rebuilt this site from Jekyll to Astro, and what I plan to write about going forward."
date: 2026-05-16
tags: ["meta", "writing"]
---

I rebuilt this site from scratch. The previous version was a hand-styled Jekyll page that had served me well since 2021, but it had become increasingly painful to add anything new — every change meant editing inline Bootstrap styles in a tangle of `<div>`s.

The new site is built with [Astro](https://astro.dev/), with the blog backed by a content collection of Markdown / MDX files. It is still hosted as static HTML on GitHub Pages, but built via GitHub Actions instead of Jekyll's old-school workflow. Math via [KaTeX](https://katex.org/) and code syntax highlighting via [Shiki](https://shiki.style/) both work out of the box, which makes writing physics and code notes much friendlier.


