# Minimal JSON parser/serializer shared by the test suite and the MCP server.
# Kept dependency-free so we don't pull JSON3 / JSON into the main project.
#
# This file is `include`d by:
#   - mcp/server.jl
#   - test/benchmark.jl
#   - test/runtests.jl (for Python-CIBSA cross-validation fixtures)
#
# The functions defined here are added to the global scope of the includer.

function skip_ws(s, i)
    while i <= length(s) && s[i] in (' ', '\t', '\n', '\r')
        i += 1
    end
    return i
end

function parse_json_value(s, i)
    i = skip_ws(s, i)
    i > length(s) && error("Unexpected end of JSON")
    c = s[i]
    if c == '"'
        return parse_json_string(s, i)
    elseif c == '{'
        return parse_json_object(s, i)
    elseif c == '['
        return parse_json_array(s, i)
    elseif c == 't'
        return (true, i + 4)
    elseif c == 'f'
        return (false, i + 5)
    elseif c == 'n'
        return (nothing, i + 4)
    elseif isdigit(c) || c == '-'
        return parse_json_number(s, i)
    else
        error("Unexpected character at position $i: '$(s[i])' ($(Int(s[i])))")
    end
end

function parse_json_string(s, i)
    i += 1  # skip opening "
    buf = IOBuffer()
    while i <= length(s) && s[i] != '"'
        if s[i] == '\\'
            i += 1
            if s[i] == 'n'; write(buf, '\n')
            elseif s[i] == 't'; write(buf, '\t')
            elseif s[i] == 'r'; write(buf, '\r')
            elseif s[i] == '"'; write(buf, '"')
            elseif s[i] == '\\'; write(buf, '\\')
            elseif s[i] == '/'; write(buf, '/')
            elseif s[i] == 'u'
                # \uXXXX Unicode escape
                hex = s[i+1:i+4]
                cp = parse(UInt32, hex; base=16)
                i += 4
                # Handle surrogate pairs
                if 0xD800 <= cp <= 0xDBFF && i + 2 <= length(s) && s[i+1] == '\\' && s[i+2] == 'u'
                    hex2 = s[i+3:i+6]
                    lo = parse(UInt32, hex2; base=16)
                    cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00)
                    i += 6
                end
                write(buf, Char(cp))
            else write(buf, s[i])
            end
        else
            write(buf, s[i])
        end
        i += 1
    end
    return (String(take!(buf)), i + 1)  # skip closing "
end

function parse_json_number(s, i)
    j = i
    if j <= length(s) && s[j] == '-'
        j += 1
    end
    while j <= length(s) && isdigit(s[j])
        j += 1
    end
    is_float = false
    if j <= length(s) && s[j] == '.'
        is_float = true
        j += 1
        while j <= length(s) && isdigit(s[j])
            j += 1
        end
    end
    if j <= length(s) && s[j] in ('e', 'E')
        is_float = true
        j += 1
        if j <= length(s) && s[j] in ('+', '-')
            j += 1
        end
        while j <= length(s) && isdigit(s[j])
            j += 1
        end
    end
    numstr = s[i:j-1]
    val = is_float ? parse(Float64, numstr) : parse(Int, numstr)
    return (val, j)
end

function parse_json_array(s, i)
    i += 1  # skip [
    arr = Any[]
    i = skip_ws(s, i)
    if i <= length(s) && s[i] == ']'
        return (arr, i + 1)
    end
    while true
        val, i = parse_json_value(s, i)
        push!(arr, val)
        i = skip_ws(s, i)
        if i > length(s) || s[i] == ']'
            return (arr, i + 1)
        end
        i += 1  # skip comma
    end
end

function parse_json_object(s, i)
    i += 1  # skip {
    obj = Dict{String, Any}()
    i = skip_ws(s, i)
    if i <= length(s) && s[i] == '}'
        return (obj, i + 1)
    end
    while true
        i = skip_ws(s, i)
        key, i = parse_json_string(s, i)
        i = skip_ws(s, i)
        i += 1  # skip colon
        val, i = parse_json_value(s, i)
        obj[key] = val
        i = skip_ws(s, i)
        if i > length(s) || s[i] == '}'
            return (obj, i + 1)
        end
        i += 1  # skip comma
    end
end

function parse_json_file(path::String)
    return parse_json_value(read(path, String), 1)[1]
end

# ── Serializer ──────────────────────────────────────────────────────────────

function json_serialize(v::AbstractString)
    buf = IOBuffer()
    write(buf, '"')
    for c in v
        if c == '"';       write(buf, "\\\"")
        elseif c == '\\';  write(buf, "\\\\")
        elseif c == '\n';  write(buf, "\\n")
        elseif c == '\r';  write(buf, "\\r")
        elseif c == '\t';  write(buf, "\\t")
        elseif codepoint(c) < 0x20 || codepoint(c) >= 0x80
            cp = codepoint(c)
            if cp <= 0xFFFF
                write(buf, "\\u$(lpad(string(cp, base=16), 4, '0'))")
            else
                cp -= 0x10000
                hi = 0xD800 + (cp >> 10)
                lo = 0xDC00 + (cp & 0x3FF)
                write(buf, "\\u$(lpad(string(hi, base=16), 4, '0'))")
                write(buf, "\\u$(lpad(string(lo, base=16), 4, '0'))")
            end
        else
            write(buf, c)
        end
    end
    write(buf, '"')
    return String(take!(buf))
end

json_serialize(v::Integer) = string(v)
json_serialize(v::AbstractFloat) = string(v)
json_serialize(v::Bool) = v ? "true" : "false"
json_serialize(::Nothing) = "null"

function json_serialize(v::AbstractVector)
    return "[" * join([json_serialize(x) for x in v], ",") * "]"
end

function json_serialize(v::Tuple)
    return "[" * join([json_serialize(x) for x in v], ",") * "]"
end

function json_serialize(v::AbstractDict)
    pairs = [json_serialize(string(k)) * ":" * json_serialize(val) for (k, val) in v]
    return "{" * join(pairs, ",") * "}"
end
