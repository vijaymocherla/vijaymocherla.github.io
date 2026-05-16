# vijaymocherla.github.io

Personal website of **Sai Vijay Mocherla** — built with [Astro](https://astro.build/)
and deployed to [GitHub Pages](https://pages.github.com/) via GitHub Actions.

Live at <https://vijaymocherla.github.io>.

## Stack

- **Astro 4** — static site generator
- **MDX** — Markdown + components for blog posts
- **KaTeX** — math rendering (server-side)
- **Shiki** — code syntax highlighting (light + dark themes)
- **GitHub Actions** — build & deploy

## Local development

Requires **Node 20+** and **npm**.

```bash
# install dependencies
npm install

# start the dev server at http://localhost:4321
npm run dev

# build the production site to ./dist
npm run build

# preview the built site locally
npm run preview
```

## Writing a blog post

1. Create a new Markdown or MDX file under `src/content/blog/`. Filename
   convention: `YYYY-MM-DD-slug.md`. The filename (without the date prefix
   is fine too — Astro uses the filename as the URL slug).
2. Add frontmatter:

   ```yaml
   ---
   title: "My new post"
   description: "Short summary used in <meta> and the post list."
   date: 2026-05-16
   tags: ["tag-one", "tag-two"]
   draft: false        # optional; set true to hide from build
   updated: 2026-05-20 # optional
   heroImage: /img/foo.webp  # optional; used for OG image and at top of post
   ---
   ```

3. Write the post body in Markdown. For math, use `$inline$` and `$$block$$`
   syntax. For embedded React-style components, rename to `.mdx`.
4. Commit and push to `master`. The GitHub Action will rebuild and deploy.

## Adding a research project

Edit `src/pages/research.astro` directly. Each project is an `<article class="project">`
block — duplicate one of the existing entries and adjust.

## Project structure

```
.
├── .github/workflows/deploy.yml   # build & deploy on push
├── public/                        # static assets served at /
│   ├── img/                       # photos, favicons, project figures
│   └── docs/                      # CV and PDF reports
├── src/
│   ├── components/                # shared UI components
│   ├── content/
│   │   ├── config.ts              # blog collection schema
│   │   └── blog/                  # blog posts (.md / .mdx)
│   ├── layouts/                   # BaseLayout, BlogPost
│   ├── pages/                     # one .astro file per route
│   │   ├── index.astro            # /
│   │   ├── research.astro         # /research/
│   │   ├── 404.astro              # /404
│   │   ├── rss.xml.js             # /rss.xml
│   │   ├── blog/[...slug].astro   # /blog/<slug>/
│   │   ├── blog/index.astro       # /blog/
│   │   ├── tags/[tag].astro       # /tags/<tag>/
│   │   └── tags/index.astro       # /tags/
│   ├── styles/global.css          # tokens (light/dark), base typography
│   └── consts.ts                  # site metadata, nav, social links
├── astro.config.mjs
├── package.json
├── tsconfig.json
└── _legacy/                       # previous Jekyll site (kept for reference)
```

## One-time GitHub Pages setup

After the first push, switch GitHub Pages to use Actions:

> **Repo → Settings → Pages → Build and deployment → Source = "GitHub Actions"**

The workflow at `.github/workflows/deploy.yml` will then build on every push
to `master` and publish to <https://vijaymocherla.github.io>.

## License

Source code: MIT. Written content: © Sai Vijay Bhaskar Mocherla.
