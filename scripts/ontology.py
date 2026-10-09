"""Parse active GO terms and propagate through is_a / part_of only."""
from functools import lru_cache

def read_go(obo):
    terms, alt, current, version = {}, {}, None, "unknown"
    def commit(term):
        if term and "id" in term and not term.get("obsolete", False):
            terms[term["id"]] = term
            for other in term["alt"]:
                alt[other] = term["id"]
    with obo.open(encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith("data-version: "):
                version = line.split(": ", 1)[1]
            if line.startswith("["):
                commit(current)
                current = {"parents": set(), "alt": []} if line == "[Term]" else None
            elif current is not None:
                if line.startswith(("id: ", "name: ", "namespace: ")):
                    k, v = line.split(": ", 1)
                    current[k] = v
                elif line.startswith("alt_id: "):
                    current["alt"].append(line.split()[1])
                elif line.startswith("is_a: "):
                    current["parents"].add(line.split()[1])
                elif line.startswith("relationship: part_of "):
                    current["parents"].add(line.split()[2])
                elif line == "is_obsolete: true":
                    current["obsolete"] = True
        commit(current)
    if not terms:
        raise ValueError("No GO terms found in ontology.")
    @lru_cache(maxsize=None)
    def ancestors(go):
        go = alt.get(go, go)
        if go not in terms:
            return frozenset()
        result = {go}
        for parent in terms[go]["parents"]:
            result.update(ancestors(alt.get(parent, parent)))
        return frozenset(result)
    return terms, alt, ancestors, version
