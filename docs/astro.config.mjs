// @ts-check
import { defineConfig } from "astro/config";
import starlight from "@astrojs/starlight";
import starlightLinksValidator from "starlight-links-validator";
import generatedPages from "./src/integrations/generated-pages.mjs";

// https://astro.build/config
export default defineConfig({
    // Project GitHub Pages of VIBTOBIlab/MethylSeq: https://vibtobilab.github.io/MethylSeq/
    site: "https://vibtobilab.github.io",
    base: "/MethylSeq/",
    integrations: [
        // Writes the parameters, changelog and citations pages before Starlight reads the content
        generatedPages(),
        starlight({
            title: "MethylSeq",
            description: "Nextflow pipeline for the analysis of bisulfite sequencing data (WGBS, RRBS and cfRRBS)",
            // TOBI lab logo (from https://vibtobilab.github.io/), with a white B on the dark theme
            logo: {
                light: "./src/assets/tobi-logo-light.png",
                dark: "./src/assets/tobi-logo-dark.png",
                alt: "TOBI lab",
            },
            favicon: "/favicon.svg",
            social: [{ icon: "github", label: "GitHub", href: "https://github.com/VIBTOBIlab/MethylSeq" }],
            customCss: ["./src/styles/custom.css"],
            // Fails the build on broken links and anchors
            plugins: [starlightLinksValidator()],
            sidebar: [
                { label: "Home", slug: "" },
                { label: "Usage", slug: "usage" },
                { label: "Parameters", slug: "parameters" },
                { label: "Output", slug: "output" },
                {
                    label: "About",
                    items: [
                        { label: "Citations", slug: "citations" },
                        { label: "Changelog", slug: "changelog" },
                        {
                            label: "Original nf-core/methylseq",
                            link: "https://nf-co.re/methylseq/2.6.0",
                            attrs: { target: "_blank" },
                        },
                    ],
                },
            ],
        }),
    ],
});
