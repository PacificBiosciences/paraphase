use crate::VString;

use std::fmt::Formatter;

#[derive(
    Default, Clone, PartialEq, Ord, PartialOrd, Eq, Hash, Copy, serde::Serialize, serde::Deserialize,
)]
pub struct VStr<'a> {
    data: &'a [u8],
}

impl<'a> VStr<'a> {
    /// # Errors
    /// `std::str::Utf8Error` if `self.data` is not valid UTF-8.
    pub fn as_str(&self) -> Result<&str, std::str::Utf8Error> {
        std::str::from_utf8(self.data)
    }
}

impl<'a> std::fmt::Display for VStr<'a> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        let view = std::str::from_utf8(self.data).map_err(|_| std::fmt::Error {})?;
        write!(f, "{view}")
    }
}

impl<'a> std::fmt::Debug for VStr<'a> {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        let view = std::str::from_utf8(self.data).map_err(|_| std::fmt::Error {})?;
        write!(f, "{view}")
    }
}

impl<'a> From<&'a VString> for VStr<'a> {
    fn from(value: &'a VString) -> Self {
        let data = &value[..];
        Self { data }
    }
}

impl<'a> PartialEq<&str> for VStr<'a> {
    fn eq(&self, other: &&str) -> bool {
        self.data == other.as_bytes()
    }
}
impl<'a> PartialOrd<&str> for VStr<'a> {
    fn partial_cmp(&self, other: &&str) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(other.as_bytes())
    }
}

impl<'a> PartialEq<String> for VStr<'a> {
    fn eq(&self, other: &String) -> bool {
        self.data == other.as_bytes()
    }
}
impl<'a> PartialOrd<String> for VStr<'a> {
    fn partial_cmp(&self, other: &String) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(other.as_bytes())
    }
}

impl<'a> PartialEq<char> for VStr<'a> {
    fn eq(&self, other: &char) -> bool {
        self.len() == 1 && self[0] as char == *other
    }
}
impl<'a> PartialOrd<char> for VStr<'a> {
    fn partial_cmp(&self, other: &char) -> Option<std::cmp::Ordering> {
        self.partial_cmp(&other.to_string())
    }
}
impl<'a> PartialEq<u8> for VStr<'a> {
    fn eq(&self, other: &u8) -> bool {
        self.len() == 1 && self[0] == *other
    }
}
impl<'a> PartialOrd<u8> for VStr<'a> {
    fn partial_cmp(&self, other: &u8) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(&[*other][..])
    }
}
impl<'a> PartialEq<[u8]> for VStr<'a> {
    fn eq(&self, other: &[u8]) -> bool {
        self.data.eq(other)
    }
}
impl<'a> PartialOrd<[u8]> for VStr<'a> {
    fn partial_cmp(&self, other: &[u8]) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(other)
    }
}
impl<'a, const N: usize> PartialEq<[u8; N]> for VStr<'a> {
    fn eq(&self, other: &[u8; N]) -> bool {
        self.data.eq(&other[..])
    }
}
impl<'a, const N: usize> PartialOrd<[u8; N]> for VStr<'a> {
    fn partial_cmp(&self, other: &[u8; N]) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(&other[..])
    }
}
impl<'a, const N: usize> PartialEq<&[u8; N]> for VStr<'a> {
    fn eq(&self, other: &&[u8; N]) -> bool {
        self.data.eq(&other[..])
    }
}
impl<'a, const N: usize> PartialOrd<&[u8; N]> for VStr<'a> {
    fn partial_cmp(&self, other: &&[u8; N]) -> Option<std::cmp::Ordering> {
        self.data.partial_cmp(&other[..])
    }
}

impl<'a> From<&VStr<'a>> for VStr<'a> {
    fn from(value: &VStr<'a>) -> Self {
        *value
    }
}

impl<'a> From<VStr<'a>> for &'a str {
    fn from(val: VStr<'a>) -> Self {
        std::str::from_utf8(&val).expect("Failed to unwrap VStr.")
    }
}

impl<'a> From<&'a String> for VStr<'a> {
    fn from(value: &'a String) -> Self {
        let data = value.as_bytes();
        Self { data }
    }
}

impl<'a> From<&'a str> for VStr<'a> {
    fn from(value: &'a str) -> Self {
        let data = value.as_bytes();
        Self { data }
    }
}

impl<'a> From<&&'a str> for VStr<'a> {
    fn from(value: &&'a str) -> Self {
        let data = value.as_bytes();
        Self { data }
    }
}
impl<'a> From<&'a [u8]> for VStr<'a> {
    fn from(value: &'a [u8]) -> Self {
        Self { data: value }
    }
}
impl<'a> From<&&'a [u8]> for VStr<'a> {
    fn from(value: &&'a [u8]) -> Self {
        Self { data: value }
    }
}
impl<'a> std::ops::Deref for VStr<'a> {
    type Target = &'a [u8];
    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        use crate::VString;
        let v = VStr::from("GATTACA");
        assert_eq!(String::from("GATTACA"), format!("{v}"));
        assert_eq!(v[0], b'G');
        let v2 = VString::from("GATTACA");
        assert_eq!(String::from("GATTACA"), format!("{v2}"));
        assert_eq!(v2[0], b'G');
        assert_eq!(v2.vstr(), v);
    }

    /// Also check creation of `VStr` from `String` and from `&str`.
    #[test]
    fn ord_ok() {
        let mut x: std::collections::BTreeMap<VStr, i32> = std::collections::BTreeMap::new();
        x.insert(VStr::from("HELLO"), 1);
        x.insert(VStr::from("GOODBYE"), 2);
        assert!("HELLO" > "GOODBYE");
        assert!("HELLO" > "GOODBYE");
        assert!("HELLO" < "gOODBYE");
        assert_eq!(VStr::from(&String::from("Hello")), VStr::from("Hello"));
        // Make sure we can compare to chars and bytes.
        assert_eq!(VStr::from("C"), 'C');
        assert_eq!(VStr::from("C"), b'C');
        assert!(VStr::from("C") != 'G');
        assert!(VStr::from("C") != 'c');
        assert!(VStr::from("cC") != 'c');
    }

    #[test]
    fn hashmap_ok() {
        let mut x = std::collections::HashMap::<VStr, i32>::new();
        x.insert(VStr::from("HELLO"), 1);
        x.insert(VStr::from("GOODBYE"), 2);
        assert!("HELLO" > "GOODBYE");
        assert_eq!(*x.get(&VStr::from("HELLO")).unwrap(), 1);
    }

    #[test]
    fn deser_ok() {
        use crate::{VStr, VString};
        let input = VString::from("John Doe");
        let data = bincode::serialize(&input.vstr()).unwrap();
        let data_vstring = bincode::serialize(&input).unwrap();
        let res: VStr<'_> = bincode::deserialize(&data).expect("Fail");
        assert_eq!(data.len(), 16);
        assert_eq!(data_vstring.len(), 16);
        assert_eq!(data, data_vstring);
        assert_eq!(res, input);
        drop(input);
        assert_eq!(res, "John Doe");
    }
}
