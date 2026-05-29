pub use crate::vstr::VStr;

use std::fmt::Formatter;

#[derive(
    Default, Clone, PartialEq, Ord, PartialOrd, Eq, Hash, serde::Serialize, serde::Deserialize,
)]
pub struct VString {
    data: Vec<u8>,
}

impl VString {
    #[must_use]
    pub fn as_bytes(&self) -> &[u8] {
        &self.data
    }
    /// # Errors
    /// `std::str::Utf8Error` if data is not valid utf-8.
    pub fn as_str(&self) -> Result<&str, std::str::Utf8Error> {
        std::str::from_utf8(&self.data)
    }
    #[must_use]
    pub fn vstr(&self) -> VStr<'_> {
        VStr::from(&self.data[..])
    }
}

impl<'a> std::iter::FromIterator<&'a u8> for VString {
    fn from_iter<I: IntoIterator<Item = &'a u8>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.push(*x);
        }
        Self::from(ret)
    }
}

impl<'a> std::iter::FromIterator<&'a str> for VString {
    fn from_iter<I: IntoIterator<Item = &'a str>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.extend_from_slice(x.as_bytes());
        }
        Self::from(ret)
    }
}

impl<'a> std::iter::FromIterator<&'a [u8]> for VString {
    fn from_iter<I: IntoIterator<Item = &'a [u8]>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.extend_from_slice(x);
        }
        Self::from(ret)
    }
}

impl<'a> std::iter::FromIterator<&'a Vec<u8>> for VString {
    fn from_iter<I: IntoIterator<Item = &'a Vec<u8>>>(iter: I) -> Self {
        iter.into_iter().map(|x| &x[..]).collect::<Self>()
    }
}

impl std::iter::FromIterator<Vec<u8>> for VString {
    fn from_iter<I: IntoIterator<Item = Vec<u8>>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.extend_from_slice(&x);
        }
        Self::from(ret)
    }
}

impl<'a> std::iter::FromIterator<&'a String> for VString {
    fn from_iter<I: IntoIterator<Item = &'a String>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.extend_from_slice(x.as_bytes());
        }
        Self::from(ret)
    }
}

impl std::iter::FromIterator<String> for VString {
    fn from_iter<I: IntoIterator<Item = String>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.extend_from_slice(x.as_bytes());
        }
        Self::from(ret)
    }
}

impl std::iter::FromIterator<u8> for VString {
    fn from_iter<I: IntoIterator<Item = u8>>(iter: I) -> Self {
        let iter = iter.into_iter();
        let mut ret: Vec<u8> = Vec::with_capacity(iter.size_hint().0);
        for x in iter {
            ret.push(x);
        }
        Self::from(ret)
    }
}

impl From<&str> for VString {
    fn from(value: &str) -> Self {
        let data = value.as_bytes().to_vec();
        Self { data }
    }
}
impl From<&[u8]> for VString {
    fn from(value: &[u8]) -> Self {
        Self::from(value.to_vec())
    }
}

impl<const N: usize> From<[u8; N]> for VString {
    fn from(value: [u8; N]) -> Self {
        Self::from(value.to_vec())
    }
}

impl<const N: usize> From<&[u8; N]> for VString {
    fn from(value: &[u8; N]) -> Self {
        Self::from(value.to_vec())
    }
}

impl From<VStr<'_>> for VString {
    fn from(value: VStr<'_>) -> Self {
        Self::from(&value[..])
    }
}
impl From<&VStr<'_>> for VString {
    fn from(value: &VStr<'_>) -> Self {
        Self::from(*value)
    }
}

impl From<&String> for VString {
    fn from(value: &String) -> Self {
        Self::from(value.bytes().collect::<Vec<_>>())
    }
}
impl From<&VString> for VString {
    fn from(value: &VString) -> Self {
        value.clone()
    }
}
impl From<String> for VString {
    fn from(value: String) -> Self {
        Self::from(value.bytes().collect::<Vec<_>>())
    }
}

impl From<Vec<u8>> for VString {
    fn from(value: Vec<u8>) -> Self {
        Self { data: value }
    }
}

impl<'a> From<&'a VString> for &'a str {
    fn from(val: &'a VString) -> Self {
        val.as_str().unwrap_or_else(|e| {
            panic!(
                "Failed to decode VString. Error: {e}. Data: {:?}",
                &val.data[..]
            )
        })
    }
}

impl std::fmt::Display for VString {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        let view = std::str::from_utf8(&self.data).map_err(|_| std::fmt::Error {})?;
        write!(f, "{view}")
    }
}

impl std::fmt::Debug for VString {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), std::fmt::Error> {
        let view = std::str::from_utf8(&self.data).map_err(|_| std::fmt::Error {})?;
        write!(f, "{view}")
    }
}

impl std::ops::Deref for VString {
    type Target = Vec<u8>;
    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl std::ops::DerefMut for VString {
    fn deref_mut(&mut self) -> &mut Vec<u8> {
        &mut self.data
    }
}

impl PartialOrd<VStr<'_>> for VString {
    fn partial_cmp(&self, other: &VStr) -> Option<std::cmp::Ordering> {
        self.vstr().partial_cmp(other)
    }
}
impl PartialOrd<VString> for VStr<'_> {
    fn partial_cmp(&self, other: &VString) -> Option<std::cmp::Ordering> {
        let vs = other.vstr();
        (*self).partial_cmp(&vs)
    }
}

impl PartialEq<VStr<'_>> for VString {
    fn eq(&self, other: &VStr) -> bool {
        self.vstr().eq(other)
    }
}
impl PartialEq<VString> for VStr<'_> {
    fn eq(&self, other: &VString) -> bool {
        *self == other.vstr()
    }
}

impl PartialEq<&str> for VString {
    fn eq(&self, other: &&str) -> bool {
        self.vstr().eq(other)
    }
}
impl PartialOrd<&str> for VString {
    fn partial_cmp(&self, other: &&str) -> Option<std::cmp::Ordering> {
        self.vstr().partial_cmp(other)
    }
}
impl PartialEq<[u8]> for VString {
    fn eq(&self, other: &[u8]) -> bool {
        self.vstr().eq(other)
    }
}
impl PartialOrd<[u8]> for VString {
    fn partial_cmp(&self, other: &[u8]) -> Option<std::cmp::Ordering> {
        self.vstr().partial_cmp(other)
    }
}
impl PartialEq<String> for VString {
    fn eq(&self, other: &String) -> bool {
        self.vstr().eq(other)
    }
}
impl PartialOrd<String> for VString {
    fn partial_cmp(&self, other: &String) -> Option<std::cmp::Ordering> {
        self.vstr().partial_cmp(other)
    }
}

impl PartialEq<char> for VString {
    fn eq(&self, other: &char) -> bool {
        self.len() == 1 && self[0] as char == *other
    }
}
impl PartialOrd<char> for VString {
    fn partial_cmp(&self, other: &char) -> Option<std::cmp::Ordering> {
        self.partial_cmp(&other.to_string())
    }
}
impl PartialEq<u8> for VString {
    fn eq(&self, other: &u8) -> bool {
        self.len() == 1 && self[0] == *other
    }
}
impl PartialOrd<u8> for VString {
    fn partial_cmp(&self, other: &u8) -> Option<std::cmp::Ordering> {
        self.data[..].partial_cmp(&[*other][..])
    }
}
impl<const N: usize> PartialEq<[u8; N]> for VString {
    fn eq(&self, other: &[u8; N]) -> bool {
        self.data.eq(&other[..])
    }
}
impl<const N: usize> PartialOrd<[u8; N]> for VString {
    fn partial_cmp(&self, other: &[u8; N]) -> Option<std::cmp::Ordering> {
        self.data[..].partial_cmp(&other[..])
    }
}
impl<const N: usize> PartialEq<&[u8; N]> for VString {
    fn eq(&self, other: &&[u8; N]) -> bool {
        self.data.eq(&other[..])
    }
}
impl<const N: usize> PartialOrd<&[u8; N]> for VString {
    fn partial_cmp(&self, other: &&[u8; N]) -> Option<std::cmp::Ordering> {
        self.data[..].partial_cmp(&other[..])
    }
}

impl std::borrow::Borrow<[u8]> for VString {
    fn borrow(&self) -> &[u8] {
        &self[..]
    }
}
impl std::borrow::Borrow<[u8]> for VStr<'_> {
    fn borrow(&self) -> &[u8] {
        &self[..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let v = VString::from("GATTACA");
        assert_eq!(String::from("GATTACA"), format!("{v}"));
        assert_eq!(v[0], b'G');
        let mut v = VString::default();
        v.push(b'A');
        v.push(b'T');
        v.push(b'C');
        assert_eq!(&format!("{v}")[..], "ATC");
    }
    #[test]
    fn constructors_ok() {
        use crate::VString;
        let seq = VString::from("GATTACA");
        let seq2 = VString::from(&b"GATTACA"[..]);
        let seq3 = VString::from(String::from("GATTACA"));
        let seq4 = VString::from(&String::from("GATTACA"));
        assert_eq!(seq, seq2);
        assert_eq!(seq, seq3);
        assert_eq!(seq, seq4);
    }
    #[test]
    fn collect_vstring_ok() {
        let seq = b"ACGT"
            .iter()
            .enumerate()
            .map(|(count, item)| (count + 1, item))
            .collect::<Vec<_>>();
        let mut tmpvec: Vec<u8> = Vec::new();
        for (count, byte) in seq {
            for _ in 0..count {
                tmpvec.push(*byte);
            }
        }
        let res = tmpvec.iter().collect::<VString>();
        assert_eq!(VString::from(&b"ACCGGGTTTT"[..]), res);
    }
    #[test]
    fn ord_ok() {
        let mut x: std::collections::BTreeMap<VString, i32> = std::collections::BTreeMap::new();
        x.insert(VString::from("HELLO"), 1);
        x.insert(VString::from("GOODBYE"), 2);
        let hello = VString::from("HELLO");
        assert!(hello > "GOODBYE");
        assert!(hello > "GOODBYE");
        assert!(hello > "GOODBYE");
        assert!(hello < "gOODBYE");
        // Make sure we can compare to chars and bytes.
        assert_eq!(VString::from("C"), 'C');
        assert_eq!(VString::from("C"), b'C');
    }
    #[test]
    fn het_eq_ok() {
        let x = VString::from("HELLO");
        let y = VStr::from("HELLO");
        assert_eq!(x, y);
        assert_eq!(y, x);
        let y = VStr::from("GOODBYE");
        assert!(x != y);
        assert!(y != x);
        assert_eq!(y, [b'G', b'O', b'O', b'D', b'B', b'Y', b'E']);
        assert_eq!(y, &[b'G', b'O', b'O', b'D', b'B', b'Y', b'E']);
        assert_eq!(y, b"GOODBYE");
    }
    #[test]
    fn hashmap_ok() {
        let mut x = std::collections::HashMap::<VString, i32>::new();
        x.insert(VString::from("HELLO"), 1);
        x.insert(VString::from("GOODBYE"), 2);
        let hello = VString::from("HELLO");
        assert!(hello > "GOODBYE");
        assert_eq!(x.get(&hello), Some(1).as_ref());
    }
    #[test]
    fn borrow_ok() {
        let mut x = std::collections::HashMap::<VString, i32>::new();
        x.insert(VString::from("HELLO"), 1);
        x.insert(VString::from("GOODBYE"), 2);
        assert_eq!(x.get(&b"HELLO"[..]), Some(&1i32));
        let mut x = std::collections::HashMap::<VStr<'_>, i32>::new();
        x.insert(VStr::from("HELLO"), 1);
        x.insert(VStr::from("GOODBYE"), 2);
        assert_eq!(x.get(&b"HELLO"[..]), Some(&1i32));
    }
}
