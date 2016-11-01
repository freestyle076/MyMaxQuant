/// 
/// Copyright (c) 2008 Max Planck Institute of Biochemistry,
/// Department of Proteomics and Signal Transduction. All rights reserved.
/// 
using System.Collections.Generic;

namespace MaxQuant.Util {
	/// <summary>
	/// Represents a set of objects (without repetition) of the same type.
	/// </summary>
	/// <typeparam name="T">The type of the objects to be stored.</typeparam>
	public class HashSet<T> {
		/// <summary>
		/// The set elements are stored as keys in this dictionary. The values are left empty.
		/// </summary>
		private readonly Dictionary<T, object> dict = new Dictionary<T, object>();

		/// <summary>
		/// Initializes the set with the members of the given list.
		/// </summary>
		/// <param name="list">List containing the elements for initialization.</param>
		public HashSet(IEnumerable<T> list) {
			AddAll(list);
		}

		/// <summary>
		/// Creates an empty set.
		/// </summary>
		public HashSet() {
		}

		/// <summary>
		/// Returns the number of elements.
		/// </summary>
		public int Count {
			get { return dict.Count; }
		}

		/// <summary>
		/// Checks whether the given object is an element of the set.
		/// </summary>
		/// <param name="x">The object to be checked whether or not it is an element.</param>
		/// <returns>Returns true iff the given object is an element of the set.</returns>
		public bool Contains(T x) {
			return dict.ContainsKey(x);
		}
		
		/// <summary>
		/// Adds the given object to the set. It does nothing if it is already contained.
		/// </summary>
		/// <param name="x">Object to be added.</param>
		public void Add(T x) {
			if (dict.ContainsKey(x)) {
				return;
			}
			dict.Add(x, null);
		}

		/// <summary>
		/// Removes the given element from the set. It does nothing if it is not contained.
		/// </summary>
		/// <param name="x">Object to be removed.</param>
		public void Remove(T x) {
			if (!dict.ContainsKey(x)) {
				return;
			}
			dict.Remove(x);
		}

		/// <summary>
		/// Returns an array containing all elements in no particular order.
		/// </summary>
		/// <returns>Array with elements.</returns>
		public T[] ToArray() {
			return DataUtil.GetKeys(dict);
		}

		/// <summary>
		/// Adds all members of the given <code>IEnumerable</code>
		/// </summary>
		/// <param name="x">List of elements to be added.</param>
		public void AddAll(IEnumerable<T> x) {
			foreach (T y in x) {
				Add(y);
			}
		}

		/// <summary>
		/// Returns an element. Returns <code>default(T)</code> if empty.
		/// </summary>
		/// <returns>An element.</returns>
		public T GetAny() {
			if(Count == 0) {
				return default(T);
			}
			return DataUtil.GetKeys(dict)[0];
		}

		/// <summary>
		/// Removes all elements.
		/// </summary>
		public void Clear() {
			dict.Clear();
		}
	}
}