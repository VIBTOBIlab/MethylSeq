/**
 * Astro integration that generates the pages built from files of the pipeline:
 *
 * - parameters.md: the template src/templates/parameters.md, with the placeholder
 *   `<!-- schema-parameters -->` replaced by one table per section of
 *   nextflow_schema.json, so the page is always in sync with the schema used to
 *   validate the parameters.
 * - changelog.md and citations.md: CHANGELOG.md and CITATIONS.md of the repository.
 *
 * The pages are written to src/content/docs/ (ignored by git) when the dev server
 * starts and before each build. The dev server restarts when a source file changes.
 *
 * Author: Edoardo Giuili (edoardo.giuili@ugent.be)
 */

import { readFileSync, writeFileSync } from "node:fs";
import { fileURLToPath } from "node:url";

const PLACEHOLDER = "<!-- schema-parameters -->";
const REPO = new URL("../../../", import.meta.url);
const DOCS = new URL("../content/docs/", import.meta.url);

const SOURCES = {
    schema: new URL("nextflow_schema.json", REPO),
    template: new URL("../templates/parameters.md", import.meta.url),
    changelog: new URL("CHANGELOG.md", REPO),
    citations: new URL("CITATIONS.md", REPO),
};

/** Escape a text of the schema, e.g. '<OUTDIR>' would be read as an HTML tag. */
const escapeHtml = (text) => String(text).replaceAll("&", "&amp;").replaceAll("<", "&lt;").replaceAll(">", "&gt;");

/** Make a text safe for a Markdown table cell. */
const cell = (text) => String(text).replaceAll("|", "\\|").replaceAll("\n", " ").trim();

/** Description of a parameter, with its help text and allowed values. */
function describe(spec) {
    const parts = [escapeHtml(spec.description ?? "")];
    if (spec.help_text) parts.push(escapeHtml(spec.help_text));
    if (spec.enum) parts.push("Allowed values: " + spec.enum.map((value) => `\`${value}\``).join(", "));
    if (spec.minimum != null) parts.push(`Minimum: \`${spec.minimum}\`.`);
    return parts.filter(Boolean).map(cell).join("<br>");
}

/** Default value of a parameter, if any. */
const defaultValue = (spec) => ("default" in spec ? `\`${spec.default}\`` : "");

/** Type of a parameter, e.g. 'boolean, string' for a parameter with several types. */
const typeOf = (spec) => [spec.type ?? ""].flat().join(", ");

/** Markdown table of the visible parameters of a schema section. */
function sectionTable(section) {
    const required = new Set(section.required ?? []);
    const rows = Object.entries(section.properties ?? {})
        .filter(([, spec]) => !spec.hidden)
        .map(
            ([name, spec]) =>
                `| \`--${name}\` | ${describe(spec)} | \`${typeOf(spec)}\` ` +
                `| ${defaultValue(spec)} | ${required.has(name) ? "yes" : ""} |`,
        );
    if (!rows.length) return "";
    return [
        "| Parameter | Description | Type | Default | Required |",
        "| --- | --- | --- | --- | :---: |",
        ...rows,
    ].join("\n");
}

/** Schema sections in the order of 'allOf', then any section missing from it. */
function sections(schema) {
    // Draft-07 schemas (nf-validation) use 'definitions', newer ones '$defs'
    const definitions = schema.definitions ?? schema.$defs ?? {};
    const ordered = (schema.allOf ?? []).map((ref) => ref.$ref.split("/").at(-1));
    ordered.push(...Object.keys(definitions).filter((name) => !ordered.includes(name)));
    return ordered.filter((name) => name in definitions).map((name) => definitions[name]);
}

/** Markdown with one heading and table per section of the schema. */
function renderSchema() {
    const schema = JSON.parse(readFileSync(SOURCES.schema, "utf8"));
    const blocks = [];
    for (const section of sections(schema)) {
        const table = sectionTable(section);
        if (!table) continue;
        blocks.push(`## ${section.title ?? ""}`);
        if (section.description) blocks.push(section.description);
        blocks.push(table);
    }
    return blocks.join("\n\n");
}

/** A Markdown file of the repository as a page, with its first heading replaced by a title. */
function repoPage(source, title, description) {
    const body = readFileSync(source, "utf8").replace(/^# .*\n+/, "");
    return `---\ntitle: ${title}\ndescription: ${description}\n---\n\n${body}`;
}

function generatePages() {
    const parameters = readFileSync(SOURCES.template, "utf8").replace(PLACEHOLDER, () => renderSchema());
    writeFileSync(new URL("parameters.md", DOCS), parameters);
    writeFileSync(
        new URL("changelog.md", DOCS),
        repoPage(SOURCES.changelog, "Changelog", "The changes of each version of the pipeline."),
    );
    writeFileSync(
        new URL("citations.md", DOCS),
        repoPage(SOURCES.citations, "Citations", "The tools used by the pipeline and their publications."),
    );
}

export default function generatedPages() {
    return {
        name: "methylseq-generated-pages",
        hooks: {
            "astro:config:setup": ({ addWatchFile }) => {
                generatePages();
                for (const source of Object.values(SOURCES)) addWatchFile(fileURLToPath(source));
            },
        },
    };
}
