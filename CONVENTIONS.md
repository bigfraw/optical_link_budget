# Documentation rules

All documentation in this repository uses ASD-STE100 Simplified Technical
English (STE). This applies to module docstrings, function docstrings, code
comments, README files, and commit messages.

## Rules

- Write short sentences. Use a maximum of 20 words for an instruction. Use a
  maximum of 25 words for a description.
- Write one instruction in one sentence.
- Use the active voice.
- Use the present tense when it is possible.
- Start an instruction with a verb.
- Use simple and approved words. Use one word for one meaning. Use one term for
  one technical concept.
- Keep the articles (a, an, the). Do not remove them.
- Do not use the -ing form of a verb, unless it is part of a technical name.
- Do not use idioms, slang, or figurative language.
- Write units, symbols, and names in a consistent form.

## Examples

Not STE:
"Left alone, that self-normalisation silently cancels the beam-broadening loss
out of the result -- easily ~10 dB, not a rounding error."

STE:
"The kernel normalises the irradiance to its own short-term waist. This removes
the beam-broadening loss. For a 600 km uplink the loss is approximately 10 dB.
The code corrects this."

Not STE:
"coupled_flux.py vendors the Dios kernels so we kill the external dep."

STE:
"olb/turbulence/coupled_flux.py holds the Dios coupled-flux kernels. olb copied
them from my_analysis_modules. olb no longer imports that repository."

## For agents

Every subagent prompt that writes code or documentation must include these
rules. Add an "ASD-STE100" instruction block to the prompt.
