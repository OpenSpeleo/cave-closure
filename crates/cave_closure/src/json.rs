//! Strict JSON parsing, including recursive duplicate-key rejection.
use crate::*;
use serde::{
    Deserialize, Deserializer,
    de::{self, MapAccess, SeqAccess, Visitor},
};
use serde_json::Value;
use std::fmt;

struct Unique(Value);
impl<'de> Deserialize<'de> for Unique {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct Strict;
        impl<'de> Visitor<'de> for Strict {
            type Value = Unique;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("JSON with unique object keys")
            }
            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Unique, A::Error> {
                let mut result = serde_json::Map::new();
                while let Some((key, value)) = map.next_entry::<String, Unique>()? {
                    if result.insert(key.clone(), value.0).is_some() {
                        return Err(de::Error::custom(format!("duplicate object key {key}")));
                    }
                }
                Ok(Unique(Value::Object(result)))
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Unique, A::Error> {
                let mut result = Vec::new();
                while let Some(v) = seq.next_element::<Unique>()? {
                    result.push(v.0);
                }
                Ok(Unique(Value::Array(result)))
            }
            fn visit_bool<E: de::Error>(self, v: bool) -> Result<Unique, E> {
                Ok(Unique(Value::Bool(v)))
            }
            fn visit_str<E: de::Error>(self, v: &str) -> Result<Unique, E> {
                Ok(Unique(Value::String(v.into())))
            }
            fn visit_string<E: de::Error>(self, v: String) -> Result<Unique, E> {
                Ok(Unique(Value::String(v)))
            }
            fn visit_i64<E: de::Error>(self, v: i64) -> Result<Unique, E> {
                Ok(Unique(v.into()))
            }
            fn visit_u64<E: de::Error>(self, v: u64) -> Result<Unique, E> {
                Ok(Unique(v.into()))
            }
            fn visit_f64<E: de::Error>(self, v: f64) -> Result<Unique, E> {
                serde_json::Number::from_f64(v)
                    .map(|v| Unique(Value::Number(v)))
                    .ok_or_else(|| de::Error::custom("non-finite number"))
            }
            fn visit_unit<E: de::Error>(self) -> Result<Unique, E> {
                Ok(Unique(Value::Null))
            }
            fn visit_none<E: de::Error>(self) -> Result<Unique, E> {
                Ok(Unique(Value::Null))
            }
        }
        deserializer.deserialize_any(Strict)
    }
}
fn reject_positional(value: &Value, schema: &Value, root: &Value) -> Result<(), Error> {
    if let Some(reference) = schema.get("$ref").and_then(Value::as_str)
        && let Some(target) = reference.strip_prefix('#').and_then(|p| root.pointer(p))
    {
        return reject_positional(value, target, root);
    }
    if schema.get("type").and_then(Value::as_str) == Some("object") && !value.is_object() {
        return Err(invalid(
            "expected JSON object, positional arrays are not accepted",
        ));
    }
    if let (Some(object), Some(properties)) = (
        value.as_object(),
        schema.get("properties").and_then(Value::as_object),
    ) {
        for (key, v) in object {
            if let Some(s) = properties.get(key) {
                reject_positional(v, s, root)?;
            }
        }
    }
    if let (Some(array), Some(items)) = (value.as_array(), schema.get("items")) {
        for v in array {
            reject_positional(v, items, root)?;
        }
    }
    for keyword in ["oneOf", "anyOf", "allOf"] {
        if let Some(variants) = schema.get(keyword).and_then(Value::as_array) {
            // Tagged unions select by the literal kind before descending.
            for variant in variants {
                let properties = variant.get("properties");
                let discriminator = properties
                    .and_then(|p| p.get("kind"))
                    .and_then(|s| s.get("const"));
                if discriminator.is_none() || discriminator == value.get("kind") {
                    reject_positional(value, variant, root)?;
                }
            }
        }
    }
    Ok(())
}
pub fn from_json(input: &str) -> Result<NormalizedProblem, Error> {
    let mut deserializer = serde_json::Deserializer::from_str(input);
    let value = Unique::deserialize(&mut deserializer)?.0;
    deserializer.end()?;
    let schema = serde_json::to_value(schemars::schema_for!(NormalizedProblem))?;
    reject_positional(&value, &schema, &schema)?;
    let problem: NormalizedProblem = serde_json::from_value(value)?;
    if problem.schema_version != 1 {
        return Err(invalid("unsupported schema version"));
    }
    Ok(problem)
}
pub fn input_schema() -> Value {
    serde_json::to_value(schemars::schema_for!(NormalizedProblem)).expect("schema serializes")
}
pub fn result_schema() -> Value {
    serde_json::to_value(schemars::schema_for!(SolveReport)).expect("schema serializes")
}
